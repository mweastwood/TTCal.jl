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
                sources::Vector{PointSource})
    # 1. Discard sources that are below the horizon.
    sources = filter(sources) do source
        isabovehorizon(frame,source)
    end

    # 2. Sort the sources in order of decreasing flux.
    #sort!(sources, by=source->source.I, rev=true)

    # 3. For each source:
    for i = 1:length(sources)
        name = sources[i].name
        l,m = getlm(frame,sources[i])
        lold = l; mold = m # TEMP

        # a) Subtract all other sources.
        other_sources = PointSource[]
        for j = 1:length(sources)
            i == j && continue
            push!(other_sources,sources[j])
        end
        subtracted = subsrc(frame,data,u,v,w,ν,other_sources)

        # b) Solve for the source's position
        l,m = fitvis_lm(subtracted,flags,l,m,u,v,w,ν,ant1,ant2)
        az,el = lm2azel(l,m)

        # c) Solve for the source's spectrum
        I,Q,U,V,reffreq,index = fitvis_spec(subtracted,flags,l,m,u,v,w,ν,ant1,ant2)

        #println("----")
        #@show name lold,l mold,m sources[i].I,I sources[i].Q,Q sources[i].U,U sources[i].V,V sources[i].index,index
        sources[i] = PointSource(name,Direction("AZEL",az,el),
                                 I,Q,U,V,reffreq,index)
    end
    sources
end

function fitvis_lm(data::Array{Complex64,3},
                   flags::Array{Bool,3},
                   l::Float64,
                   m::Float64,
                   u::Vector{quantity(Float64,Meter)},
                   v::Vector{quantity(Float64,Meter)},
                   w::Vector{quantity(Float64,Meter)},
                   ν::Vector{quantity(Float64,Hertz)},
                   ant1::Vector{Int32},
                   ant2::Vector{Int32})
    Nfreq = length(ν)
    Nbase = length(u)
    λ = [c/ν[β] for β = 1:Nfreq]

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
    for β = 1:Nfreq
        F[β] /= count[β]
        dF[:,β] /= count[β]
        d2F[:,:,β] /= count[β]
    end

    # Calculate how far to step in each direction
    dl = 0.0
    dm = 0.0
    normalization = 0.0
    for β = 1:Nfreq
        count[β] .== 0 && continue
        δ = slice(d2F,:,:,β) \ slice(dF,:,β)
        weight = count[β]/(λ[β]/Meter)
        dl += δ[1] * weight
        dm += δ[2] * weight
        normalization += weight
    end
    dl /= normalization
    dm /= normalization
    force_to_horizon(l-dl,m-dm)
end

function fitvis_spec(data::Array{Complex64,3},
                     flags::Array{Bool,3},
                     l::Float64,
                     m::Float64,
                     u::Vector{quantity(Float64,Meter)},
                     v::Vector{quantity(Float64,Meter)},
                     w::Vector{quantity(Float64,Meter)},
                     ν::Vector{quantity(Float64,Hertz)},
                     ant1::Vector{Int32},
                     ant2::Vector{Int32})
    I,Q,U,V = getspec_internal(data,flags,l,m,u,v,w,ν,ant1,ant2)

    # power law fit
    mask = (!isnan(I)) .* (I .> 0)
    reffreq = 47e6Hertz
    if sum(mask) > 0.5length(mask)
        # Give a reasonable anser if there are too many channels masked
        α = [0.0 0.0]
    else
        x = log10(ν[mask]/reffreq)
        y = log10(I[mask])
        α = [ones(length(x)) x]\y
    end

    # polarization fractions
    # (this is messing things up somehow)
    fQ = mean(Q[mask]./I[mask])
    fU = mean(U[mask]./I[mask])
    fV = mean(V[mask]./I[mask])

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

