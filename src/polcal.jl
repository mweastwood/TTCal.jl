# Copyright (c) 2015 Michael Eastwood
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

################################################################################
# Public Interface

function polcal(ms::Table,
                sources::Vector{PointSource},
                criteria::StoppingCriteria;
                force_imaging_columns::Bool = false,
                model_already_present::Bool = false)
    frame = reference_frame(ms)
    dir   = phase_dir(ms)
    u,v,w = uvw(ms)
    ν = freq(ms)
    ant1,ant2 = ants(ms)

    sources = filter(source -> isabovehorizon(frame,source),sources)

    Nfreq = length(ν)
    Nant  = numrows(Table(ms[kw"ANTENNA"]))

    gains = ones(Complex128,2,2,Nant,Nfreq)
    gain_flags = zeros(Bool,Nant,Nfreq)

    data  = Tables.checkColumnExists(ms,"CORRECTED_DATA")? ms["CORRECTED_DATA"] : ms["DATA"]
    model = model_already_present? ms["MODEL_DATA"] : genvis(dir,sources,u,v,w,ν)
    data_flags = ms["FLAG"]
    row_flags  = ms["FLAG_ROW"]
    for α = 1:length(row_flags)
        if row_flags[α]
            data_flags[:,:,α] = true
        end
    end

    if force_imaging_columns || Tables.checkColumnExists(ms,"MODEL_DATA")
        ms["MODEL_DATA"] = model
    end

    polcal!(gains,gain_flags,data,model,data_flags,
            ant1,ant2,criteria)
    gains,gain_flags
end

################################################################################
# Internal Interface

function polcal!(gains::Array{Complex128,4},
                 gain_flags::Array{Bool,2},
                 data::Array{Complex64,3},
                 model::Array{Complex64,3},
                 data_flags::Array{Bool,3},
                 ant1::Vector{Int32},
                 ant2::Vector{Int32},
                 criteria::StoppingCriteria)
    # Transpose the data and model arrays to create a better memory access pattern
    data  = permutedims(data, (1,3,2))
    model = permutedims(model,(1,3,2))
    data_flags = permutedims(data_flags,(1,3,2))

    Nfreq = size(gains,4)
    for β = 1:Nfreq
        polcal_onechannel!(slice(gains,:,:,:,β),
                           slice(gain_flags,:,β),
                           slice( data,:,:,β),
                           slice(model,:,:,β),
                           slice(data_flags,:,:,β),
                           ant1, ant2,
                           criteria)
    end
end

function polcal_onechannel!(gains, gain_flags,
                            data, model, data_flags,
                            ant1::Vector{Int32},
                            ant2::Vector{Int32},
                            criteria::StoppingCriteria)
    Nant = size(gains,3)
    Nbase = size(data,2)

    # If the entire channel is flagged, don't bother calibrating.
    # Just flag the solutions and move on.
    if all(data_flags)
        gain_flags[:] = true
        return
    end

    for α = 1:Nbase
        if any(data_flags[:,α])
            data[:,α] = 0
            model[:,α] = 0
        end
    end

    # This is an especially good approximation if we've
    # already applied a bandpass calibration.
    for i = 1:Nant
        gains[:,:,i] = [1 0;
                        0 1]
    end

    oldχ2 = typemax(Float64)
    iter = 0
    converged = false
    while !converged && iter < criteria.maxiter
        χ2,gain_flags[:] = newtonstep!(gains,data,model,1/100)
        if abs(χ2 - oldχ2)/abs(oldχ2) < criteria.tolerance
            converged = true
        end
        oldχ2 = χ2
        iter += 1
    end

    nothing
end

function newtonstep!(gains,data,model,damping_factor)
    Nant = size(gains,3)
    Nbase = div(Nant*(Nant-1),2)
    # Calculate χ² and its derivatives
    χ   = 0.0
    dχ  = zeros(8Nant)
    d2χ = zeros(8Nant) # We will approximate d2χ (the Hessian) as diagonal
    for ant1 = 1:Nant, ant2 = ant1+1:Nant
        α = div((ant1-1)*(2Nant-ant1+2),2) + ant2 - ant1 + 1
        V = [data[1,α] data[2,α];
             data[3,α] data[4,α]]
        M = [model[1,α] model[2,α];
             model[3,α] model[4,α]]
        G1 = slice(gains,:,:,ant1)
        G2 = slice(gains,:,:,ant2)

        # Compute the residual
        res = reinterpret(Float64,vec(V - G1*M*G2'))

        # Derivative of the model
        D = zeros(Complex128,2,2)
        dM1 = zeros(Float64,8,8)
        dM2 = zeros(Float64,8,8)
        @inbounds for p = 0:7
            # (with respect to antenna 1 parameters)
            D[:] = 0
            D[div(p,2)+1] = ifelse(mod(p,2)==0,1.0+0.0im,0.0+1.0im)
            X = D*M*G2'
            dM1[1,p+1] = real(X[1,1])
            dM1[2,p+1] = imag(X[1,1])
            dM1[3,p+1] = real(X[2,1])
            dM1[4,p+1] = imag(X[2,1])
            dM1[5,p+1] = real(X[1,2])
            dM1[6,p+1] = imag(X[1,2])
            dM1[7,p+1] = real(X[2,2])
            dM1[8,p+1] = imag(X[2,2])

            # (with respect to antenna 2 parameters)
            D[:] = 0
            D[div(p,2)+1] = ifelse(mod(p,2)==0,1.0+0.0im,0.0+1.0im)
            X = G1*M*D'
            dM2[1,p+1] = real(X[1,1])
            dM2[2,p+1] = imag(X[1,1])
            dM2[3,p+1] = real(X[2,1])
            dM2[4,p+1] = imag(X[2,1])
            dM2[5,p+1] = real(X[1,2])
            dM2[6,p+1] = imag(X[1,2])
            dM2[7,p+1] = real(X[2,2])
            dM2[8,p+1] = imag(X[2,2])
        end

        # Calculate the corresponding derivatives of χ²
        χ += sum(abs2(res))
        @inbounds for p1 = 0:7
            for p2 = 0:7
                dχ[8ant1-7+p1] -= res[p2+1].*dM1[p2+1,p1+1]
                dχ[8ant2-7+p1] -= res[p2+1].*dM2[p2+1,p1+1]
                # The following two lines give the only nonzero terms
                # along the diagonal of the Hessian
                d2χ[8ant1-7+p1] += dM1[p2+1,p1+1].*dM1[p2+1,p1+1]
                d2χ[8ant2-7+p1] += dM2[p2+1,p1+1].*dM2[p2+1,p1+1]
            end
        end
    end
    step = -(diagm(d2χ)+χ*damping_factor*I)\dχ
    flags = isnan(step)
    step[flags] = 0
    step_reinterpreted = reshape(reinterpret(Complex128,step),(2,2,Nant))
    gains[:] = gains + step_reinterpreted
    χ,flags[1:8:end]
end

