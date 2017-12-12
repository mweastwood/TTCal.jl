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

struct Source{Shapes<:AbstractVector{<:AbstractShape}}
    name   :: String
    shapes :: Shapes
end

function Source(name::String, shape::AbstractShape)
    Source(name, [shape])
end

function isabovehorizon(frame::ReferenceFrame, source::Source; threshold=0)
    for shape in source.shapes
        if !isabovehorizon(frame, shape, threshold=threshold)
            return false
        end
    end
    true
end

function isrising(frame::ReferenceFrame, source::Source)
    for shape in source.shapes
        if !isrising(frame, shape)
            return false
        end
    end
    true
end

function total_flux(source, ν)
    output = zero(StokesVector)
    for shape in source.shapes
        output += shape.spectrum(ν)
    end
    output
end

function mean_direction(frame, source, ν)
    output = Measures.UnnormalizedDirection(dir"ITRF", 0, 0, 0)
    for shape in source.shapes
        stokes  = shape.spectrum(ν)
        itrf    = measure(frame, shape.direction, dir"ITRF")
        output += stokes.I*itrf
    end
    Direction(output) # normalizes
end

function Base.:*(number::Real, source::Source)
    Source(source.name, [number*shape for shape in source.shapes])
end

Base.:*(source::Source, number::Real) = number*source

