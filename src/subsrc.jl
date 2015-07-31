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

################################################################################
# Public Interface

function subsrc!(ms::Table,sources::Vector{PointSource})
    phase_dir = MeasurementSets.phase_direction(ms)
    u,v,w = MeasurementSets.uvw(ms)
    ν = MeasurementSets.frequency(ms)

    frame = ReferenceFrame()
    set!(frame,MeasurementSets.position(ms))
    set!(frame,MeasurementSets.time(ms))
    sources = filter(source -> isabovehorizon(frame,source),sources)

    data  = MeasurementSets.corrected_data(ms)

    subtracted = subsrc(frame,phase_dir,data,u,v,w,ν,sources)
    ms["CORRECTED_DATA"] = subtracted
    subtracted
end

@doc """
Subtract all of the measured flux from a given direction.
This should be used to remove RFI, preventing it from contaminating
the other routines.
""" ->
function subsrc!(ms::Table,dir::Direction)
    phase_dir = MeasurementSets.phase_direction(ms)
    u,v,w = MeasurementSets.uvw(ms)
    ν = MeasurementSets.frequency(ms)
    ant1,ant2 = MeasurementSets.antennas(ms)

    frame = ReferenceFrame()
    set!(frame,MeasurementSets.position(ms))
    set!(frame,MeasurementSets.time(ms))

    data  = MeasurementSets.corrected_data(ms)
    flags = MeasurementSets.flags(ms)

    j2000 = measure(frame,dir,Measures.J2000)
    l,m = dir2lm(frame,phase_dir,dir)
    I,Q,U,V = getspec_internal(data,flags,l,m,u,v,w,ν,ant1,ant2)

    model = zeros(Complex64,4,length(ν),length(u))
    genvis!(model,I,Q,U,V,l,m,u,v,w,ν)
    subtracted = model # (just a rename, subtracted and model are the same array)
    for i in eachindex(model)
        subtracted[i] = data[i] - model[i]
    end
    ms["CORRECTED_DATA"] = subtracted
    subtracted
end

################################################################################
# Internal Interface

function subsrc(frame::ReferenceFrame,
                phase_dir::Direction,
                data::Array{Complex64,3},
                u::Vector{Float64},
                v::Vector{Float64},
                w::Vector{Float64},
                ν::Vector{Float64},
                sources::Vector{PointSource})
    sources = filter(sources) do source
        isabovehorizon(frame,source)
    end

    model = genvis(frame,phase_dir,sources,u,v,w,ν)
    # Re-use the space allocated for the model visibilities
    # to store the model subtracted visibilities.
    for i in eachindex(model)
        model[i] = data[i]-model[i]
    end
    model
end

