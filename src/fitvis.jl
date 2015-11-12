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
    fitvis(ms::MeasurementSet, direction::Direction;
           maxiter = 20, tolerance = 1e-3, minuvw = 0.0) -> l,m

Fit for the location of a point source near the given direction.
"""
function fitvis(ms::MeasurementSet,
                direction::Direction;
                maxiter::Int = 20,
                tolerance::Float64 = 1e-3,
                minuvw::Float64 = 0.0)
    if !isabovehorizon(ms.frame,direction)
        error("Direction is below the horizon.")
    end

    j2000 = measure(ms.frame,direction,dir"J2000")
    l,m   = direction_cosines(ms.phase_direction,j2000)

    data  = get_corrected_data(ms)
    flags = get_flags(ms)
    flag_short_baselines!(flags,minuvw,ms.u,ms.v,ms.w,ms.ν)

    fitvis(data,flags,l,m,
           ms.u,ms.v,ms.w,ms.ν,
           ms.ant1,ms.ant2,
           maxiter,tolerance)
end

function fitvis(data,flags,l,m,
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
    F   = 0.0        # flux
    dF  = [0.0, 0.0] # [ ∂F/∂l ∂F/∂m ]
    d2F = [0.0 0.0;  # [ ∂²F/∂l²  ∂²F/∂l∂m
           0.0 0.0]  #   ∂²F/∂l∂m ∂²F/∂m² ]
    count = 0 # the number of baselines used in the calculation

    πi = π*1im
    λ = [c/ν[β] for β = 1:Nfreq]
    @inbounds for α = 1:Nbase
        ant1[α] == ant2[α] && continue # don't use auto-correlations
        for β = 1:Nfreq
            any(slice(flags,:,β,α)) && continue
            # use the Stokes I visibilities only
            V = 0.5*(data[1,β,α]+data[4,β,α])*conj(fringe[β,α])
            uλ = u[α]/λ[β]
            vλ = v[α]/λ[β]
            F        += real(V)
            dF[1]    += real(V * -2πi * uλ)
            dF[2]    += real(V * -2πi * vλ)
            d2F[1,1] += real(V * (-2πi)^2 * uλ^2)
            d2F[2,1] += real(V * (-2πi)^2 * uλ*vλ)
            d2F[1,2] += real(V * (-2πi)^2 * uλ*vλ)
            d2F[2,2] += real(V * (-2πi)^2 * vλ^2)
            count += 1
        end
    end

    # Normalizing isn't strictly necessary because we care
    # about the ratio of the Hessian to the gradient, but
    # we do it here just to avoid confusion later in life
    F /= count
    dF /= count
    d2F /= count

    # Calculate how far to step in each direction
    δ = -d2F \ dF
    dl = δ[1]
    dm = δ[2]
    l′,m′ = force_to_horizon(l+dl,m+dm)

    # If the step size is larger than a resolution element
    # we will stay put because the fitting process is about
    # to diverge.
    n  = sqrt(1-l^2-m^2)
    n′ = sqrt(1-l′^2-m′^2)
    dθ = acos(l*l′+m*m′+n*n′)
    baseline_length = [sqrt(u[α]^2+v[α]^2+w[α]^2) for α = 1:Nbase]
    resolution = minimum(λ) / maximum(baseline_length)
    if dθ > resolution
        l′ = l
        m′ = m
    end

    [l′-l, m′-m]
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

