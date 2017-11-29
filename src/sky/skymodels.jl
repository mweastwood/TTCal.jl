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

abstract type AbstractSource end

#function Base.show(io::IO, source::RFISource)
#    print(io, "RFISource(\"", source.name, "\", ", source.position, "\", ", source.spectrum, ")")
#end
#
#function ==(lhs::PointSource, rhs::PointSource)
#    lhs.name == rhs.name && lhs.direction == rhs.direction && lhs.spectrum == rhs.spectrum
#end
#
#function ==(lhs::GaussianSource, rhs::GaussianSource)
#    lhs.name == rhs.name && lhs.direction == rhs.direction && lhs.spectrum == rhs.spectrum &&
#        lhs.major_fwhm == rhs.major_fwhm && lhs.minor_fwhm == rhs.minor_fwhm &&
#        lhs.position_angle == rhs.position_angle
#end
#
#function ==(lhs::MultiSource, rhs::MultiSource)
#    lhs.name == rhs.name && lhs.components == rhs.components
#end
#
#function isabovehorizon(frame::ReferenceFrame, source::RFISource, threshold = 0)
#    true
#end
#
#function isabovehorizon(frame::ReferenceFrame, source::MultiSource, threshold = 0)
#    for component in source.components
#        if !isabovehorizon(frame, component, threshold)
#            return false
#        end
#    end
#    true
#end
#
#function abovehorizon{T<:Source}(frame::ReferenceFrame, sources::Vector{T}, threshold = 0)
#    filter(sources) do source
#        isabovehorizon(frame, source, threshold)
#    end
#end
#
