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
    subsrc!(ms::MeasurementSet,
            sources::Vector{PointSource},
            beam::BeamModel)

Remove the list of sources from the measurement set.
"""
function subsrc!(ms::MeasurementSet,
                 sources::Vector{PointSource},
                 beam::BeamModel)
    sources = filter(source -> isabovehorizon(ms.frame,source),sources)
    data  = get_corrected_data(ms)
    model = genvis(ms,sources,beam)
    subsrc!(data,model)
    set_corrected_data!(ms,data)
    data
end

"""
    subsrc!(ms::MeasurementSet, dir::Direction)

Subtract all of the measured flux from a given direction.

This can be used to remove RFI sources provided they have
a known direction.
"""
function subsrc!(ms::MeasurementSet, dir::Direction)
    data  = get_corrected_data(ms)
    flags = get_flags(ms)

    j2000 = measure(ms.frame,dir,dir"J2000")
    l,m = dir2lm(ms.frame,ms.phase_dir,dir)
    xx,xy,yx,yy = getspec(data,flags,l,m,u,v,w,ν,ant1,ant2)

    model = genvis(xx,xy,yx,yy,l,m,ms.u,ms.v,ms.w,ms.ν)
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

