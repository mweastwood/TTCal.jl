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
                reference_antenna::Int = 1)
    frame = reference_frame(ms)
    u,v,w = uvw(ms)
    ν = freq(ms)
    ant1,ant2 = ants(ms)

    sources = filter(source -> isabovehorizon(frame,source),sources)

    Nfreq = length(ν)
    Nant  = numrows(Table(ms[kw"ANTENNA"]))

    gains = ones(Complex64,2,2,Nant,Nfreq)
    gain_flags = zeros(Bool,Nant,Nfreq)

    data  = Tables.checkColumnExists(ms,"CORRECTED_DATA")? ms["CORRECTED_DATA"] : ms["DATA"]
    model = genvis(frame,sources,u,v,w,ν)
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
            ant1,ant2,criteria,reference_antenna)
    gains,gain_flags
end

################################################################################
# Internal Interface

function polcal!(gains::Array{Complex64,4},
                 gain_flags::Array{Bool,2},
                 data::Array{Complex64,3},
                 model::Array{Complex64,3},
                 data_flags::Array{Bool,3},
                 ant1::Vector{Int32},
                 ant2::Vector{Int32},
                 criteria::StoppingCriteria,
                 reference_antenna::Int)
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
                           criteria,
                           reference_antenna)
    end
end

function polcal_onechannel!(gains, gain_flags,
                            data, model, data_flags,
                            ant1::Vector{Int32},
                            ant2::Vector{Int32},
                            criteria::StoppingCriteria,
                            reference_antenna::Int)
    # 1. If the entire channel is flagged, don't bother calibrating.
    #    Just flag the solutions and move on.
    if all(data_flags)
        gain_flags[:] = true
        return
    end

    # 2. Pack the visibilities into square, Hermitian matrices.
    square_data  = makesquare_polarized(data, ant1,ant2)
    square_model = makesquare_polarized(model,ant1,ant2)
    square_flags = makesquare_polarized(data_flags,ant1,ant2)

    # 3. Flag the auto-correlations.
    for i = 1:size(square_flags,1)
        square_flags[i,i] = true
    end

    # 4. Start the gains off at something sensible.
    square_data  = square_data  .* !square_flags
    square_model = square_model .* !square_flags
    # This is an especially good approximation if we've
    # already applied a bandpass calibration.
    Nant = size(gains,3)
    for i = 1:Nant
        gains[:,:,i] = [1 0;
                        0 1]
    end

    # 5. Iteratively improve that estimate.

    # Create rkstep! workspace variables.
    # Because "gains" is not a vector, it isn't serviced by the
    # current definition of RKWorkspace
    #workspace = RKWorkspace(gains,4)
    x′ = similar(gains)
    k  = [similar(gains) for i = 1:4]
    oldgains = similar(gains)

    iter = 0
    converged = false
    while !converged && iter < criteria.maxiter
        oldgains[:] = gains
        rkstep!(gains,x′,k,polcal_step!,RK4,square_data,square_model)
        if vecnorm(gains-oldgains)/vecnorm(oldgains) < criteria.tolerance
            converged = true
        end
        iter += 1
    end

    # 6. Fix the phase of the reference antenna.
    # TODO

    # 7. Flag the antennas with no unflagged data.
    bad_pols = squeeze(all(square_flags,1),1)
    bad_ants = zeros(Bool,Nant)
    for ant = 1:Nant
        if bad_pols[2ant-1] && bad_pols[2ant]
            bad_ants[ant] = true
        end
    end
    gain_flags[bad_ants] = true

    nothing
end

function makesquare_polarized(input,ant1::Vector{Int32},ant2::Vector{Int32})
    N = size(input,2)
    M = 2round(Integer,div(sqrt(1+8N)-1,2))
    output = zeros(eltype(input),M,M)
    for α = 1:N
        if ant1[α] == ant2[α]
            output[2ant1[α]-1,2ant1[α]-1] = input[1,α] # xx
            output[2ant1[α]-0,2ant1[α]-1] = input[2,α] # xy
            output[2ant1[α]-1,2ant1[α]-0] = input[3,α] # yx
            output[2ant1[α]-0,2ant1[α]-0] = input[4,α] # yy
        else
            output[2ant2[α]-1,2ant1[α]-1] = input[1,α] # x₁x₂
            output[2ant2[α]-0,2ant1[α]-1] = input[2,α] # x₁y₂
            output[2ant2[α]-1,2ant1[α]-0] = input[3,α] # y₁x₂
            output[2ant2[α]-0,2ant1[α]-0] = input[4,α] # y₁y₂

            output[2ant1[α]-1,2ant2[α]-1] = conj(input[1,α]) # x₂x₁
            output[2ant1[α]-1,2ant2[α]-0] = conj(input[2,α]) # y₂x₁
            output[2ant1[α]-0,2ant2[α]-1] = conj(input[3,α]) # x₂y₁
            output[2ant1[α]-0,2ant2[α]-0] = conj(input[4,α]) # y₂y₁
        end
    end
    output
end

function polcal_inner!(output,input,V,M)
    N = size(input,3)
    Z = similar(V)
    @inbounds for j = 1:N, i = 1:N
        # TODO: make this cleaner
        # The matrix multiplication is manually inlined to prevent the construction of
        # a temporary matrix. A cleaner way to express this would be nice.
        #Z[2i-1:2i,2j-1:2j] = sub(input,:,:,i)' * sub(M,2i-1:2i,2j-1:2j)
        Z[2i-1,2j-1] = conj(input[1,1,i]) * M[2i-1,2j-1] + conj(input[2,1,i]) * M[2i-0,2j-1]
        Z[2i-0,2j-1] = conj(input[1,2,i]) * M[2i-1,2j-1] + conj(input[2,2,i]) * M[2i-0,2j-1]
        Z[2i-1,2j-0] = conj(input[1,1,i]) * M[2i-1,2j-0] + conj(input[2,1,i]) * M[2i-0,2j-0]
        Z[2i-0,2j-0] = conj(input[1,2,i]) * M[2i-1,2j-0] + conj(input[2,2,i]) * M[2i-0,2j-0]
    end
    @inbounds for j = 1:N
        output[:,:,j] = sub(Z,:,2j-1:2j) \ sub(V,:,2j-1:2j) - sub(input,:,:,j)
    end
end
const polcal_step! = RKInnerStep{:polcal_inner!}()

