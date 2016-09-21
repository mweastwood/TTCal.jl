# Copyright (c) 2015, 2016 Michael Eastwood
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
    fitvis(visibilities, metadata, beam, source; tolerance = 1e-3)
    fitvis(visibilities, metadata, beam, direction; tolerance = 1e-3)

*Description*

Fit for the location of a source near the given direction.

*Arguments*

* `visibilities` - the visibilities measured by the interferometer
* `metadata` - the metadata describing the interferometer
* `source` - the source model whose direction we want to measure
* `direction` - alternatively we can specify a direction, in which case
    we will measure the direction to a flat spectrum point source

*Keyword Arguments*

* `tolerance` - the absolute tolerance used to test for convergence
    (defaults to `1e-3` where the units are of the direction cosine)

Note that there is currently no way to set the maximum number of iterations.
Tests with NLopt showed that if the routine exceeded the maximum number of
function evaluations, NLopt would simply return the starting value instead
of its current best-guess. This is not the behavior I wanted so I've
removed the ability to set the maximum number of iterations.
"""
function fitvis(visibilities::Visibilities, meta::Metadata, source::Source;
                tolerance::Float64 = 1e-3)
    fitvis_internal(visibilities, meta, source, tolerance)
end

function fitvis(visibilities::Visibilities, meta::Metadata, direction::Direction;
                tolerance::Float64 = 1e-3)
    flat = PowerLaw(1, 0, 0, 0, 10e6, [0.0])
    point = PointSource("dummy", direction, flat)
    fitvis(visibilities, meta, point, tolerance=tolerance)
end

function fitvis_internal(visibilities, meta, source::Source, tolerance)
    frame = reference_frame(meta)
    direction = get_mean_direction(frame, source)
    if isabovehorizon(frame, source)
        data, flags = rotate_visibilities(visibilities, meta, source)
        direction = fitvis_internal(data, flags, meta, direction, tolerance)
    end
    measure(frame, direction, dir"J2000")
end

function fitvis_internal(data, flags, meta, direction::Direction, tolerance)
    uvw = UVW(meta)

    count = 0
    start = [0.0, 0.0, 0.0]
    objective(x, g) = (count += 1; fitvis_nlopt_objective(x, g, data, flags, meta, direction, uvw))
    constraint(x, g) = fitvis_nlopt_constraint(x, g, direction)

    opt = Opt(:LN_COBYLA, 3)
    equality_constraint!(opt, constraint)
    max_objective!(opt, objective)
    initial_step!(opt, sind(10/60))
    xtol_abs!(opt, tolerance)
    flux, offset, _ = optimize(opt, start)

    Direction(dir"ITRF", direction.x + offset[1], direction.y + offset[2], direction.z + offset[3])
end

function fitvis_nlopt_objective(vector, gradient, data, flags, meta, phase_direction, uvw)
    # We have a custom visibility generation routine because we only want Stokes I visibilities.
    # Hopefully this functionality will be mostly folded into `getspec` one day.
    l = phase_direction.x + vector[1]
    m = phase_direction.y + vector[2]
    n = phase_direction.z + vector[3]
    source_direction = Direction(dir"ITRF", l, m, n)
    delays = geometric_delays(meta.antennas, source_direction, phase_direction)

    flux = 0.0
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
            flux += real(V)
            count += 1
        end
    end
    flux /= count
    flux
end

function fitvis_nlopt_constraint(vector, gradient, phase_direction)
    l = phase_direction.x + vector[1]
    m = phase_direction.y + vector[2]
    n = phase_direction.z + vector[3]
    output = l^2 + m^2 + n^2 - 1
    output
end

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

