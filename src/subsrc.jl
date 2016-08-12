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
    subsrc!(visibilities, metadata, beam, sources)
    subsrc!(visibilities, model)

*Description*

Subtract sources from the visibilities.

*Arguments*

* `visibilities` - the visibilities from which the sources will be subtracted
* `metadata` - the metadata describing the interferometer
* `beam` - the primary beam model
* `sources` - the list of sources to subtract
"""
function subsrc!{T<:Source}(visibilities::Visibilities, meta::Metadata, beam::BeamModel, sources::Vector{T})
    frame = reference_frame(meta)
    sources = abovehorizon(frame, sources)
    model = genvis(meta, sources)
    subsrc!(visibilities, model)
end

function subsrc!(visibilities::Visibilities, meta::Metadata, beam::BeamModel, source::Source)
    subsrc!(visibilities, meta, beam, [source])
end

function subsrc!(visibilities::Visibilities, model::Visibilities)
    for i in eachindex(visibilities.data, model.data)
        visibilities.data[i] -= model.data[i]
    end
    visibilities
end

"""
    putsrc!(visibilities, metadata, beam, sources)
    putsrc!(visibilities, model)

*Description*

Add sources to the visibilities.

*Arguments*

* `visibilities` - the visibilities to which the sources will be added
* `metadata` - the metadata describing the interferometer
* `beam` - the primary beam model
* `sources` - the list of sources to add
"""
function putsrc!{T<:Source}(visibilities::Visibilities, meta::Metadata, beam::BeamModel, sources::Vector{T})
    frame = reference_frame(meta)
    sources = abovehorizon(frame, sources)
    model = genvis(meta, sources)
    putsrc!(visibilities, model)
end

function putsrc!(visibilities::Visibilities, meta::Metadata, beam::BeamModel, source::Source)
    putsrc!(visibilities, meta, beam, [source])
end

function putsrc!(visibilities, model)
    for i in eachindex(visibilities.data, model.data)
        visibilities.data[i] += model.data[i]
    end
    visibilities
end

