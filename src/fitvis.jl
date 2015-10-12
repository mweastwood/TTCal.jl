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

"""
    fitvis(ms::MeasurementSet,
           sources::Vector{PointSource};
           maxiter::Int = 20,
           tolerance::Float64 = 1e-3,
           minuvw::Float64 = 0.0)

Fit for the location of each point source.
"""
function fitvis(ms::MeasurementSet,
                sources::Vector{PointSource};
                maxiter::Int = 20,
                tolerance::Float64 = 1e-3,
                minuvw::Float64 = 0.0)
    sources = filter(source -> isabovehorizon(ms.frame,source),sources)
    Nsource = length(sources)
    data  = get_corrected_data(ms)
    flags = get_flags(ms)

    # Flag all of the short baselines
    for α = 1:ms.Nbase, β = 1:ms.Nfreq
        if sqrt(ms.u[α]^2 + ms.v[α]^2 + ms.w[α]^2) < minuvw*c/ms.ν[β]
            flags[:,β,α] = true
        end
    end

    l = zeros(Nsource)
    m = zeros(Nsource)
    for i = 1:Nsource
        l′,m′ = lm(ms.frame,ms.phase_direction,sources[i])
        l[i],m[i] = fitvis_onesource(data,flags,l′,m′,
                                     ms.u,ms.v,ms.w,ms.ν,
                                     ms.ant1,ms.ant2,
                                     maxiter,tolerance)
    end
    l,m
end

function fitvis_onesource(data,flags,l,m,
                          u,v,w,ν,ant1,ant2,
                          maxiter,tolerance)
    lm = [l;m]
    converged = @iterate(FitVisStep(),RK4,maxiter,tolerance,
                         lm,data,flags,u,v,w,ν,ant1,ant2)
    lm[1],lm[2]
end

function fitvis_step(lm,data,flags,
                     u,v,w,ν,ant1,ant2)
    l = lm[1]
    m = lm[2]
    Nfreq = length(ν)
    Nbase = length(u)

    # Calculate the flux in the given direction
    # and its derivatives with respect to direction.
    fringe = fringepattern(l,m,u,v,w,ν)
    F   = zeros(Nfreq)     # flux at each frequency
    dF  = zeros(2,Nfreq)   # [ ∂F/∂l ∂F/∂m ]
    d2F = zeros(2,2,Nfreq) # [ ∂²F/∂l²  ∂²F/∂l∂m
                           #   ∂²F/∂l∂m ∂²F/∂m² ]
    count = zeros(Int,Nfreq) # The number of baselines used in the calculation

    πi = π*1im
    λ = [c/ν[β] for β = 1:Nfreq]
    @inbounds for α = 1:Nbase
        ant1[α] == ant2[α] && continue # Don't use auto-correlations
        for β = 1:Nfreq
            any(slice(flags,:,β,α)) && continue
            z = conj(fringe[β,α])
            V = 0.5*(data[1,β,α]+data[4,β,α])*z # Stokes I only
            uu = u[α]/λ[β]
            vv = v[α]/λ[β]
            F[β]       += real(V)
            dF[1,β]    += real(V * -2πi * uu)
            dF[2,β]    += real(V * -2πi * vv)
            d2F[1,1,β] += real(V * (-2πi)^2 * uu^2)
            d2F[2,1,β] += real(V * (-2πi)^2 * uu*vv)
            d2F[1,2,β] += real(V * (-2πi)^2 * uu*vv)
            d2F[2,2,β] += real(V * (-2πi)^2 * vv^2)
            count[β] += 1
        end
    end
    @inbounds for β = 1:Nfreq
        F[β] /= count[β]
        dF[:,β] /= count[β]
        d2F[:,:,β] /= count[β]
    end

    # Calculate how far to step in each direction
    # TODO: weight using `count` and wavelength appropriately
    # (higher frequencies should be given more weight)
    dl = 0.0
    dm = 0.0
    normalization = 0.0
    for β = 1:Nfreq
        count[β] .== 0 && continue
        δ = -slice(d2F,:,:,β) \ slice(dF,:,β)
        dl += δ[1]
        dm += δ[2]
        normalization += 1
    end
    dl /= normalization
    dm /= normalization

    # Don't let the source lie beyond the horizon
    l_horizon,m_horizon = force_to_horizon(l+dl,m+dm)
    [l_horizon-l;m_horizon-m]
end

immutable FitVisStep <: StepFunction end
call(::FitVisStep,lm,data,flags,u,v,w,ν,ant1,ant2) = fitvis_step(lm,data,flags,u,v,w,ν,ant1,ant2)

doc"""
    force_to_horizon(l,m)

This function forces the coordinates $(l,m)$ to be above the horizon.

Although this is a nasty hack, it is necessary for fitting
some sources that are near the horizon.
"""
function force_to_horizon(l,m)
    r2 = l^2+m^2
    if r2 > 1
        θ = atan2(l,m)
        l = (1.0-1e-12)*sin(θ)
        m = (1.0-1e-12)*cos(θ)
    end
    l,m
end

