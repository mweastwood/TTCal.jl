# Copyright (c) 2015-2017 Michael Eastwood
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

#"""
#    fitvis(visibilities, metadata, beam, source; tolerance = 1e-3)
#    fitvis(visibilities, metadata, beam, direction; tolerance = 1e-3)
#
#*Description*
#
#Fit for the location of a source near the given direction.
#
#*Arguments*
#
#* `visibilities` - the visibilities measured by the interferometer
#* `metadata` - the metadata describing the interferometer
#* `source` - the source model whose direction we want to measure
#* `direction` - alternatively we can specify a direction, in which case
#    we will measure the direction to a flat spectrum point source
#
#*Keyword Arguments*
#
#* `tolerance` - the absolute tolerance used to test for convergence
#    (defaults to `1e-3` where the units are of the direction cosine)
#
#Note that there is currently no way to set the maximum number of iterations.
#Tests with NLopt showed that if the routine exceeded the maximum number of
#function evaluations, NLopt would simply return the starting value instead
#of its current best-guess. This is not the behavior I wanted so I've
#removed the ability to set the maximum number of iterations.
#"""
function fitvis(dataset::Dataset, source::Source; tolerance::Float64 = 1e-3)
    fitvis_internal(deepcopy(dataset), source, tolerance)
end

#function fitvis(visibilities::Visibilities, meta::Metadata, direction::Direction;
#                tolerance::Float64 = 1e-3)
#    flat = PowerLaw(1, 0, 0, 0, 10e6, [0.0])
#    point = PointSource("dummy", direction, flat)
#    fitvis(visibilities, meta, point, tolerance=tolerance)
#end

function fitvis_internal(dataset, source::Source, tolerance)
    # TODO: handle multiple time integrations correctly
    frame = ReferenceFrame(dataset.metadata)
    if isabovehorizon(frame, source)
        rotate_phase_center!(dataset, source)
        direction = measure(frame, dataset.metadata.phase_centers[1], dir"ITRF")
        newdirection = fitvis_internal(frame, dataset, direction, tolerance)
        source = fitvis_rotate_source_model(frame, source, direction, newdirection)
    end
    source
end

function fitvis_internal(frame, dataset, direction::Direction, tolerance)
    count = 0
    function objective(x, g)
        count += 1
        newdirection = Direction(dir"ITRF", x[1]*u"rad", x[2]*u"rad")
        fitvis_nlopt_objective(dataset, newdirection)
    end

    opt = Opt(:LN_SBPLX, 2)
    max_objective!(opt, objective)
    long = ustrip(longitude(direction))
    lat  = ustrip( latitude(direction))
    lower_bounds!(opt, [long-deg2rad(1), lat-deg2rad(1)])
    upper_bounds!(opt, [long+deg2rad(1), lat+deg2rad(1)])
    ftol_abs!(opt, tolerance)
    start = [long, lat]
    flux, finish, _ = optimize(opt, start)

    Direction(dir"ITRF", finish[1]*u"rad", finish[2]*u"rad")
end

function fitvis_nlopt_objective(dataset, direction::Direction)
    stokes = getflux(dataset, direction)
    stokes.I
end

function fitvis_rotate_source_model(frame, source, from, to)
    R = Measures.RotationMatrix(from, to)
    shapes = [fitvis_rotate_shape(frame, shape, R) for shape in source.shapes]
    Source(source.name, shapes)
end

function fitvis_rotate_shape(frame, shape::Point, R)
    direction = R*measure(frame, shape.direction, R.sys)
    direction = measure(frame, direction, dir"J2000")
    Point(direction, shape.spectrum)
end

function fitvis_rotate_shape(frame, shape::Gaussian, R)
    direction = R*measure(frame, shape.direction, R.sys)
    direction = measure(frame, direction, dir"J2000")
    Gaussian(direction, shape.spectrum, shape.major_fwhm, shape.minor_fwhm, shape.position_angle)
end

function fitvis_rotate_shape(frame, shape::Disk, R)
    direction = R*measure(frame, shape.direction, R.sys)
    direction = measure(frame, direction, dir"J2000")
    Disk(direction, shape.spectrum, shape.radius)
end

function fitvis_rotate_shape(frame, shape::Shapelet, R)
    direction = R*measure(frame, shape.direction, R.sys)
    direction = measure(frame, direction, dir"J2000")
    Shapelet(direction, shape.spectrum, shape.scale, shape.coeff)
end

