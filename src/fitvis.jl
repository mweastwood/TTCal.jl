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

@doc """
Fit the visibilities to a model of point sources. The input model needs
to have the positions of the point sources relatively close, but the flux
can be wildly off.
""" ->
function fitvis(ms::Table,
                sources::Vector{Source},
                criteria::StoppingCriteria)
    frame = reference_frame(ms)
    data  = ms["CORRECTED_DATA"]
    flags = ms["FLAG"]
    u,v,w = uvw(ms)
    ν = freq(ms)
    ant1,ant2 = ants(ms)

    fitvis(frame,data,flags,u,v,w,ν,ant1,ant2,sources,criteria)
end

################################################################################
# Internal Interface

@doc """
1. Sort the sources in order of decreasing flux.
2. For each source, subtract all other sources.
3. Fit for the source.
""" ->
function fitvis(frame::ReferenceFrame,
                data::Array{Complex64,3},
                flags::Array{Bool,3},
                u::Vector{quantity(Float64,Meter)},
                v::Vector{quantity(Float64,Meter)},
                w::Vector{quantity(Float64,Meter)},
                ν::Vector{quantity(Float64,Hertz)},
                ant1::Vector{Int32},
                ant2::Vector{Int32},
                sources::Vector{Source},
                criteria::StoppingCriteria)
    # Flag short baselines?
    for α = 1:length(u)
        if sqrt(u[α]^2+v[α]^2+w[α]^2) < 50Meter
            flags[:,:,α] = true
        end
    end

    # 1. Discard sources that are below the horizon.
    sources = filter(sources) do source
        isabovehorizon(frame,source)
    end
    @show sources

    # 2. Sort the sources in order of decreasing flux.
    sort!(sources, by=source->source.flux, rev=true)
    @show sources

    # 3. Convert all the sources to `RFISources`
    # (this removes the spectral smoothness constraint)
    output_sources = [RFISource(frame,source,ν) for source in sources]
    @show output_sources

    # 3. For each source:
    for i = 1:length(output_sources)
        @show i
        # a) Subtract all other sources.
        other_sources = RFISource[]
        for j = 1:length(output_sources)
            i == j && continue
            push!(other_sources,output_sources[j])
        end
        subtracted = subsrc(frame,data,u,v,w,ν,other_sources)

        # b) Obtain the source's spectrum.
        l,m = lm(frame,output_sources[i])
        @time spectrum = getspec_internal(subtracted,flags,l,m,u,v,w,ν,ant1,ant2)
        xx = squeeze(spectrum[1,:],1)
        yy = squeeze(spectrum[4,:],1)
        I = 0.5*(xx + yy)

        # c) Fit for the source's position.
        @time l,m = fitvis_outer(subtracted,flags,l,m,u,v,w,ν,ant1,ant2,spectrum,criteria)
        az,el = lm2azel(l,m)
        azel  = Direction("AZEL",az,el)

        # d) Update the source.
        output_sources[i] = RFISource(output_sources[i].name,azel,I,ν)
    end
    output_sources
end

@doc """
Take Runge-Kutta steps until the position converges.
""" ->
function fitvis_outer(data::Array{Complex64,3},
                      flags::Array{Bool,3},
                      l::Float64,
                      m::Float64,
                      u::Vector{quantity(Float64,Meter)},
                      v::Vector{quantity(Float64,Meter)},
                      w::Vector{quantity(Float64,Meter)},
                      ν::Vector{quantity(Float64,Hertz)},
                      ant1::Vector{Int32},
                      ant2::Vector{Int32},
                      spectrum::Array{Float64,2},
                      criteria::StoppingCriteria)
    position = [l,m]
    oldposition = copy(position)
    workspace = RKWorkspace(position,2)

    # Take Runge-Kutta steps until the position converges.
    iter = 0
    converged = false
    while !converged && iter < criteria.maxiter
        oldposition[:] = position
        @time rkstep!(position,workspace,fitvis_step!,data,flags,u,v,w,ν,ant1,ant2,spectrum)

        @show (position[1] - oldposition[1]) (position[2] - oldposition[2])
        if abs(position[1] - oldposition[1]) < criteria.tolerance &&
           abs(position[2] - oldposition[2]) < criteria.tolerance
            converged = true
        end
        iter += 1
    end

    l = position[1]
    m = position[2]
    l,m
end

@doc """
1. Fit for the position and flux of the source at each frequency.
2. Get a single position for the source.
""" ->
function fitvis_inner!(output::Vector{Float64},
                       input::Vector{Float64},
                       data::Array{Complex64,3},
                       flags::Array{Bool,3},
                       u::Vector{quantity(Float64,Meter)},
                       v::Vector{quantity(Float64,Meter)},
                       w::Vector{quantity(Float64,Meter)},
                       ν::Vector{quantity(Float64,Hertz)},
                       ant1::Vector{Int32},
                       ant2::Vector{Int32},
                       spectrum::Array{Float64,2})
    l = input[1]
    m = input[2]
    n = sqrt(1-l^2-m^2)
    Nbase = length(u)
    Nfreq = length(ν)

    fringe = fringepattern(l,m,u,v,w,ν)
    delta  = zeros(Complex64,4Nbase)
    # Columns contain: [∂V/∂l ∂V/∂m]
    # Where V is the contribution to the visibility of a single point source,
    #       and l and m specify the position of the point source.
    matrix = Array(Complex64,4Nbase,2)
    δl = Array(Float64,Nfreq)
    δm = Array(Float64,Nfreq)

    # 1. Fit for the position of the source at each frequency.
    @inbounds for β = 1:Nfreq
        λ = c/ν[β]
        count = 0 # The number of unflagged visibilities
        for α = 1:length(u)
            ant1[α] == ant2[α] && continue
            (flags[1,β,α] || flags[2,β,α] || flags[3,β,α] || flags[4,β,α]) && continue

            # Note: constraining the complex conjugates as well forces l and m
            # to be real-valued.

            model_xx = 0.5*spectrum[1,β]*fringe[β,α]
            model_yy = 0.5*spectrum[4,β]*fringe[β,α]

            dfringe_l = 2.0im*π*(u[α].val-w[α].val*l/n)/λ.val
            dfringe_m = 2.0im*π*(v[α].val-w[α].val*m/n)/λ.val

            delta[4count+1] = data[1,β,α] - model_xx
            delta[4count+2] = data[4,β,α] - model_yy
            delta[4count+3] = conj(delta[4count+1])
            delta[4count+4] = conj(delta[4count+2])

            matrix[4count+1,1] = dfringe_l * model_xx
            matrix[4count+2,1] = dfringe_l * model_yy
            matrix[4count+3,1] = conj(matrix[4count+1,1])
            matrix[4count+4,1] = conj(matrix[4count+2,1])

            matrix[4count+1,2] = dfringe_m * model_xx
            matrix[4count+2,2] = dfringe_m * model_yy
            matrix[4count+3,2] = conj(matrix[4count+1,2])
            matrix[4count+4,2] = conj(matrix[4count+2,2])

            count += 1
        end
        params = sub(matrix,1:4count,:)\sub(delta,1:4count)
        δl[β] = real(params[1])
        δm[β] = real(params[2])
    end

    # 2. Get a single position for the source.
    # (Weight the positions by ν^2 because resolution ~ 1/ν)
    weights = ν.^2
    ΔL = sum(weights .* δl) ./ sum(weights)
    ΔM = sum(weights .* δm) ./ sum(weights)
    @show ΔL ΔM

    output[1] = ΔL
    output[2] = ΔM
    nothing
end

const fitvis_step! = RKInnerStep{:fitvis_inner!}()

@doc """
This function forces l and m to be above the horizon.
Although this is a nasty hack, it is necessary for fitting
some sources that are near the horizon.
""" ->
function force_to_horizon(l,m)
    # (check to make sure we're still above the horizon)
    r2 = l^2+m^2
    if r2 > 1
        # (move l and m back inside the horizon)
        θ = atan2(l,m)
        l = (1.0-1e-12)*sin(θ)
        m = (1.0-1e-12)*cos(θ)
    end
    l,m
end

