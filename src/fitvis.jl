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
    fitvis(visibilities, metadata, beam, source; maxiter = 20, tolerance = 1e-3)
    fitvis(visibilities, metadata, beam, direction; maxiter = 20, tolerance = 1e-3)

*Description*

Fit for the location of a point source near the given direction.

*Arguments*

*Keyword Arguments*

"""
function fitvis(visibilities::Visibilities, meta::Metadata, source::Source;
                maxiter::Int = 20, tolerance::Float64 = 1e-3)
    fitvis_internal(visibilities, meta, source, maxiter, tolerance)
end

function fitvis(visibilities::Visibilities, meta::Metadata, direction::Direction;
                maxiter::Int = 20, tolerance::Float64 = 1e-3)
    flat = PowerLaw(1, 0, 0, 0, 10e6, [0.0])
    point = PointSource("dummy", direction, flat)
    fitvis(visibilities, meta, point, maxiter=maxiter, tolerance=tolerance)
end

# Design justification
#
# The visibility measured by a given baseline is influenced by:
#     1. the flux of the source
#     2. the primary beam
#     3. the baseline fringe pattern, and
#     4. the resolved structure of the source.
# The first item is constant with respect to the source position. For the LWA
# the primary beam varies on ten degree scales whereas the refraction
# through the ionosphere is expected to occur on arcminute scales. Similarly
# the resolved structure of the source only changes as the projection of the
# baseline onto the sky changes. This will therefore also vary on ten degree
# scales.
#
# As a result of the above physical arguments we will assume only the derivative
# of the fringe pattern with respect to direction is nonzero.

function fitvis_internal(visibilities, meta, source::Source, maxiter, tolerance)
    frame = reference_frame(meta)
    direction = get_mean_direction(frame, source)
    if isabovehorizon(frame, source)
        data, flags = rotate_visibilities(visibilities, meta, source)
        direction = fitvis_internal(data, flags, meta, direction, maxiter, tolerance)
    end
    measure(frame, direction, dir"J2000")
end

function fitvis_internal(data, flags, meta, direction::Direction, maxiter, tolerance)
    vector = [0.0, 0.0, 0.0, 1.0]
    uvw = UVW(meta)
    converged = iterate(fitvisstep, RK4, maxiter, tolerance, false,
                        vector, data, flags, meta, direction, uvw)
    Direction(dir"ITRF", direction.x + vector[1], direction.y + vector[2], direction.z + vector[3])
end

function fitvis_step(vector, data, flags, meta, phase_direction, uvw)
    # We have a custom visibility generation routine because
    #     1. we want the gradient and hessian as well, and
    #     2. we only want Stokes I.
    # Hopefully this functionality will be mostly folded into `genvis` one day.
    l = phase_direction.x + vector[1]
    m = phase_direction.y + vector[2]
    n = phase_direction.z + vector[3]
    lagrange = vector[4] # the Lagrange multiplier
    source_direction = Direction(dir"ITRF", l, m, n)
    delays = geometric_delays(meta.antennas, source_direction, phase_direction)

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
            flags[α,β] && continue
            fringe = conj(fringes[antenna1]) * fringes[antenna2]
            V = data[α,β]*fringe
            uλ = uvw.u[α]/λ
            vλ = uvw.v[α]/λ
            wλ = uvw.w[α]/λ
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

    # add the contribution of the Lagrange multiplier
    gradient[1]  += lagrange * 2l
    gradient[2]  += lagrange * 2m
    gradient[3]  += lagrange * 2n
    gradient[4]   = l^2 + m^2 + n^2 - 1
    hessian[1,1] += lagrange * 2
    hessian[2,2] += lagrange * 2
    hessian[3,3] += lagrange * 2
    hessian[4,1]  = 2l
    hessian[4,2]  = 2m
    hessian[4,3]  = 2n

    # the Hessian should be symmetric
    hessian[1,2] = hessian[2,1]
    hessian[1,3] = hessian[3,1]
    hessian[1,4] = hessian[4,1]
    hessian[2,3] = hessian[3,2]
    hessian[2,4] = hessian[4,2]
    hessian[3,4] = hessian[4,3]

    # compute the step
    δ = -hessian\gradient
end

immutable FitvisStep <: StepFunction end
const fitvisstep = FitvisStep()
function call(::FitvisStep, vector, data, flags, meta, phase_direction, uvw)
    fitvis_step(vector, data, flags, meta, phase_direction, uvw)
end
return_type(::FitvisStep, vector) = Vector{Float64}

function get_mean_direction(frame, source::Source)
    measure(frame, source.direction, dir"ITRF")
end

function get_mean_direction(frame, source::MultiSource)
    x, y, z = 0.0, 0.0, 0.0
    for component in source.components
        direction = get_mean_direction(frame, component)
        x += direction.x
        y += direction.y
        z += direction.z
    end
    norm = sqrt(x^2 + y^2 + z^2)
    Direction(dir"ITRF", x/norm, y/norm, z/norm)
end

function stokes_I_only(input)
    output = zeros(Complex128, size(input))
    for idx in eachindex(input, output)
        J = input[idx]
        output[idx] = 0.5*(J.xx + J.yy)
    end
    output
end

function rotate_phase_center!(data, model)
    for idx in eachindex(data, model)
        data[idx] = data[idx] / model[idx]
    end
end

function rotate_visibilities(visibilities, meta, source)
    model = genvis_internal(meta, ConstantBeam(), [source])
    flatten_spectrum!(meta, model, source)
    stokes_I_data  = stokes_I_only(visibilities.data)
    stokes_I_model = stokes_I_only(model)
    rotate_phase_center!(stokes_I_data, stokes_I_model)
    stokes_I_data, visibilities.flags
end

