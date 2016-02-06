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
    FitvisVariables

This type is simply a container for variables that need to be calculated
once and shared between functions while fitting source positions.
"""
type FitvisVariables
    phase_center :: Direction # should always be ITRF
    u :: Vector{Float64}
    v :: Vector{Float64}
    w :: Vector{Float64}
    function FitvisVariables(meta::Metadata)
        frame = reference_frame(meta)
        phase_center = measure(frame, meta.phase_center, dir"ITRF")
        u = zeros(Nbase(meta))
        v = zeros(Nbase(meta))
        w = zeros(Nbase(meta))
        for α = 1:Nbase(meta)
            antenna1 = meta.baselines[α].antenna1
            antenna2 = meta.baselines[α].antenna2
            r1 = meta.antennas[antenna1].position
            r2 = meta.antennas[antenna2].position
            u[α] = r1.x - r2.x
            v[α] = r1.y - r2.y
            w[α] = r1.z - r2.z
        end
        new(phase_center, u, v, w)
    end
end

"""
    fitvis(visibilities::Visibilities, meta::Metadata, direction::Direction;
           maxiter = 20, tolerance = 1e-3)

Fit for the location of a point source near the given direction.
"""
function fitvis(visibilities::Visibilities, meta::Metadata, direction::Direction;
                maxiter::Int = 20, tolerance::Float64 = 1e-3)
    frame = reference_frame(meta)
    isabovehorizon(frame, direction) || error("Direction is below the horizon.")
    direction = measure(frame, direction, dir"ITRF")
    variables = FitvisVariables(meta)
    newdirection = fitvis_internal(visibilities, meta, variables, direction, maxiter, tolerance)
    measure(frame, newdirection, dir"J2000")
end

function fitvis_internal(visibilities, meta, variables, direction, maxiter, tolerance)
    vector = [direction.x, direction.y, direction.z, 1.0]
    converged = iterate(FitvisStep(), RK4, maxiter, tolerance, false,
                        vector, visibilities, meta, variables)
    Direction(dir"ITRF", vector[1], vector[2], vector[3])
end

function fitvis_step(vector, visibilities, meta, variables)
    x = vector[1]
    y = vector[2]
    z = vector[3]
    lagrange = vector[4] # the Lagrange multiplier
    direction = Direction(dir"ITRF", x, y, z)
    delays = geometric_delays(meta.antennas, direction, variables.phase_center)

    πi = π*1im
    gradient = [0.0, 0.0, 0.0, 0.0]
    hessian  = [0.0  0.0  0.0  0.0
                0.0  0.0  0.0  0.0
                0.0  0.0  0.0  0.0
                0.0  0.0  0.0  0.0]

    count = 1
    for β = 1:Nfreq(meta)
        ν = meta.channels[β]
        λ = c/ν
        fringes = delays_to_fringes(delays, ν)
        for α = 1:Nbase(meta)
            antenna1 = meta.baselines[α].antenna1
            antenna2 = meta.baselines[α].antenna2
            antenna1 == antenna2 && continue # don't use auto-correlations
            visibilities.flags[α,β] && continue
            fringe = conj(fringes[antenna1]) * fringes[antenna2]
            # use the Stokes I flux only
            V = 0.5*(visibilities.data[α,β].xx+visibilities.data[α,β].yy) * fringe
            uλ = variables.u[α]/λ
            vλ = variables.v[α]/λ
            wλ = variables.w[α]/λ
            gradient[1]  += real(V * -2πi * uλ)
            gradient[2]  += real(V * -2πi * vλ)
            gradient[3]  += real(V * -2πi * wλ)
            hessian[1,1] += real(V * (-2πi)^2 * uλ^2)
            hessian[2,1] += real(V * (-2πi)^2 * uλ*vλ)
            hessian[3,1] += real(V * (-2πi)^2 * uλ*wλ)
            hessian[2,2] += real(V * (-2πi)^2 * vλ^2)
            hessian[3,2] += real(V * (-2πi)^2 * vλ*wλ)
            hessian[3,3] += real(V * (-2πi)^2 * wλ^2)
            count += 1
        end
    end
    gradient /= count
    hessian  /= count

    # Add the contribution of the Lagrange multiplier
    gradient[1]  += lagrange * 2x
    gradient[2]  += lagrange * 2y
    gradient[3]  += lagrange * 2z
    gradient[4]   = x^2 + y^2 + z^2 - 1
    hessian[1,1] += lagrange * 2
    hessian[2,2] += lagrange * 2
    hessian[3,3] += lagrange * 2
    hessian[4,1]  = 2x
    hessian[4,2]  = 2y
    hessian[4,3]  = 2z

    # The Hessian should be symmetric
    hessian[1,2] = hessian[2,1]
    hessian[1,3] = hessian[3,1]
    hessian[1,4] = hessian[4,1]
    hessian[2,3] = hessian[3,2]
    hessian[2,4] = hessian[4,2]
    hessian[3,4] = hessian[4,3]

    # Compute the step
    δ = -hessian\gradient
end

immutable FitvisStep <: StepFunction end
call(::FitvisStep, vector, visibilities, meta, variables) = fitvis_step(vector, visibilities, meta, variables)

