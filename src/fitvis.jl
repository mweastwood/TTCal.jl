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
                sources::Vector{PointSource})
    frame = reference_frame(ms)
    data  = Tables.checkColumnExists(ms,"CORRECTED_DATA")? ms["CORRECTED_DATA"] : ms["DATA"]
    flags = ms["FLAG"]
    u,v,w = uvw(ms)
    ν = freq(ms)
    ant1,ant2 = ants(ms)

    fitvis(frame,data,flags,u,v,w,ν,ant1,ant2,sources)
end

################################################################################
# Internal Interface

function fitvis(frame::ReferenceFrame,
                data::Array{Complex64,3},
                flags::Array{Bool,3},
                u::Vector{Float64},
                v::Vector{Float64},
                w::Vector{Float64},
                ν::Vector{Float64},
                ant1::Vector{Int32},
                ant2::Vector{Int32},
                sources::Vector{PointSource})
    # 1. Discard sources that are below the horizon.
    sources = filter(sources) do source
        isabovehorizon(frame,source)
    end

    # 2. Get the flux of each source.
    # (This first pass is necessary to make sure a poor initial estimate of the source
    # flux doesn't impact the results. This can happen if the source is low in the beam
    # and is much fainter than it otherwise would be.)
    for (i,source) in enumerate(sources)
        l,m = lm(frame,source)
        I,Q,U,V,reffreq,index = fitvis_spec(data,flags,l,m,u,v,w,ν,ant1,ant2)
        sources[i] = PointSource(source.name,source.dir,I,Q,U,V,reffreq,index)
    end

    # 3. Sort by order of decreasing flux.
    # (This ensures we fit and remove the bright sources first, preventing their sidelobes
    # from contaminating the faint sources.)
    sort!(sources,by=source->source.I,rev=true)

    # 3. Fit for the position and flux of each source, subtracting each source from the
    #    data in turn.
    for (i,source) in enumerate(sources)
        l,m = getlm(frame,source)

        # a) Solve for the source's position
        l,m = fitvis_lm(data,flags,l,m,u,v,w,ν,ant1,ant2)
        dir = Direction("AZEL",lm2azel(l,m)...)

        # b) Solve for the source's spectrum
        I,Q,U,V,reffreq,index = fitvis_spec(data,flags,l,m,u,v,w,ν,ant1,ant2)

        # c) Subtract the source from the data
        sources[i] = PointSource(source.name,dir,I,Q,U,V,reffreq,index)
        data = subsrc(frame,data,u,v,w,ν,[sources[i]])
    end
    sources
end

function fitvis_lm(data::Array{Complex64,3},
                   flags::Array{Bool,3},
                   l::Float64,
                   m::Float64,
                   u::Vector{Float64},
                   v::Vector{Float64},
                   w::Vector{Float64},
                   ν::Vector{Float64},
                   ant1::Vector{Int32},
                   ant2::Vector{Int32})
    Nfreq = length(ν)
    Nbase = length(u)
    λ = [c/ν[β] for β = 1:Nfreq]

    # Calculate the angular resolution of the interferometer
    resolution = minimum(λ) / maximum(sqrt(u.^2+v.^2))

    # Calculate the flux in the given direction
    # (and its derivatives with respect to direction)
    fringe = fringepattern(l,m,u,v,w,ν)
    F   = zeros(Nfreq)     # flux at each frequency
    dF  = zeros(2,Nfreq)   # [ ∂F/∂l ∂F/∂m ]
    d2F = zeros(2,2,Nfreq) # [ ∂²F/∂l²  ∂²F/∂l∂m
                           #   ∂²F/∂l∂m ∂²F/∂m² ]
    count  = zeros(Int,Nfreq) # The number of baselines used in the calculation
    @inbounds for α = 1:Nbase
        # Don't use auto-correlations
        ant1[α] == ant2[α] && continue
        for β = 1:Nfreq
            any(slice(flags,:,β,α)) && continue
            z = conj(fringe[β,α])
            V = 0.5*(data[1,β,α]+data[4,β,α])*z
            uu = u[α]/λ[β]
            vv = v[α]/λ[β]
            F[β]       += real(V)
            dF[1,β]    += real(V * -2im*π * uu)
            dF[2,β]    += real(V * -2im*π * vv)
            d2F[1,1,β] += real(V * (-2im*π)^2 * uu^2)
            d2F[2,1,β] += real(V * (-2im*π)^2 * uu*vv)
            d2F[1,2,β] += real(V * (-2im*π)^2 * uu*vv)
            d2F[2,2,β] += real(V * (-2im*π)^2 * vv^2)
            count[β] += 1
        end
    end
    @inbounds for β = 1:Nfreq
        F[β] /= count[β]
        dF[:,β] /= count[β]
        d2F[:,:,β] /= count[β]
    end

    # Calculate how far to step in each direction
    # (weight by 1/wavelength because resolution scales as 1/wavelength)
    dl = 0.0
    dm = 0.0
    normalization = 0.0
    for β = 1:Nfreq
        count[β] .== 0 && continue
        if det(slice(d2F,:,:,β)) > 0
            # If the determinant of the Hessian is positive,
            # take a Newton step towards where the gradient is zero.
            δ = -slice(d2F,:,:,β) \ slice(dF,:,β)
        else
            # If the determinant of the Hessian is negative,
            # the above step will take us away from the flux maximum.
            # Therefore we'll just step in the direction of the
            # gradient.
            δ = slice(dF,:,β) * resolution
        end
        weight = count[β]/(λ[β]/Meter)
        dl += δ[1] * weight
        dm += δ[2] * weight
        normalization += weight
    end
    dl /= normalization
    dm /= normalization

    # Restrict the step size to one resolution element
    jump_size  = sqrt(dl^2+dm^2)
    if jump_size > resolution
        correction = resolution/jump_size
        dl *= correction
        dm *= correction
    end

    # Furthermore, don't let the source lie beyond the horizon
    force_to_horizon(l+dl,m+dm)
end

function fitvis_spec(data::Array{Complex64,3},
                     flags::Array{Bool,3},
                     l::Float64,
                     m::Float64,
                     u::Vector{Float64},
                     v::Vector{Float64},
                     w::Vector{Float64},
                     ν::Vector{Float64},
                     ant1::Vector{Int32},
                     ant2::Vector{Int32})
    I,Q,U,V = getspec_internal(data,flags,l,m,u,v,w,ν,ant1,ant2)

    # power law fit
    mask = (!isnan(I)) .* (I .> 0)
    reffreq = 47e6
    if sum(mask) < 5
        # Give a reasonable answer if there are too many channels masked
        α = [0.0;0.0]
        fQ = 0.0
        fU = 0.0
        fV = 0.0
    else
        x = log10(ν[mask]/reffreq)
        y = log10(I[mask])
        α = [ones(length(x)) x]\y
        # polarization fractions
        fQ = mean(Q[mask]./I[mask])
        fU = mean(U[mask]./I[mask])
        fV = mean(V[mask]./I[mask])
    end

    flux = 10^α[1]
    flux,fQ*flux,fU*flux,fV*flux,reffreq,α[2:end]
end

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

