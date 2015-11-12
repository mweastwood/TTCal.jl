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
    subsrc!(ms::MeasurementSet, sources::Vector{Source}, beam::BeamModel)

Remove the list of sources from the measurement set.
"""
function subsrc!(ms::MeasurementSet, sources::Vector{Source}, beam::BeamModel)
    sources = abovehorizon(ms.frame,sources)
    data  = get_corrected_data(ms)
    model = genvis(ms,sources,beam)
    subsrc!(data,model)
    set_corrected_data!(ms,data)
    data
end

"""
    subsrc!(ms::MeasurementSet, dir::Direction; minuvw = 0.0)

Subtract all of the measured flux from a given direction.

This can be used to remove RFI sources provided they have
a known direction.
"""
function subsrc!(ms::MeasurementSet, dir::Direction;
                 minuvw::Float64 = 0.0)
    data  = get_corrected_data(ms)
    flags = get_flags(ms)

    flag_short_baselines!(flags,minuvw,ms.u,ms.v,ms.w,ms.ν)

    j2000 = measure(ms.frame,dir,dir"J2000")
    l,m   = direction_cosines(ms.phase_direction,j2000)
    flux  = getspec(data,flags,l,m,ms.u,ms.v,ms.w,ms.ν,ms.ant1,ms.ant2)

    model = genvis(flux,l,m,ms.u,ms.v,ms.w,ms.ν)
    subsrc!(data,model)
    set_corrected_data!(ms,data)
    data
end

function subsrc!(data,model)
    for i in eachindex(data,model)
        data[i] -= model[i]
    end
    data
end

function putsrc!(data,model)
    for i in eachindex(data,model)
        data[i] += model[i]
    end
    data
end

