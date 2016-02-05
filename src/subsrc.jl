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
    subsrc!(visibilities::Visibilities, meta::Metadata, sources::Vector{Source})

Remove the list of sources from the measurement set.
"""
function subsrc!(visibilities::Visibilities, meta::Metadata, sources::Vector{Source})
    frame = reference_frame(meta)
    sources = abovehorizon(frame, sources)
    model = genvis(meta, sources)
    subsrc!(visibilities, model)
end

function subsrc!(visibilities, model)
    for i in eachindex(visibilities.data, model.data)
        visibilities.data[i] -= model.data[i]
    end
    visibilities
end

function putsrc!(visibilities, model)
    for i in eachindex(visibilities.data, model.data)
        visibilities.data[i] += model.data[i]
    end
    visibilities
end

