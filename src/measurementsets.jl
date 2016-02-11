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

immutable Antenna
    position :: Position
end

immutable Baseline
    antenna1 :: Int
    antenna2 :: Int
end

type Metadata
    antennas  :: Vector{Antenna}
    baselines :: Vector{Baseline}
    channels  :: Vector{Float64}
    phase_center :: Direction
    time :: Epoch
    beam :: BeamModel
end

"""
    collect_metadata(ms::Table, beam)

Read the interferometer's instrumental parameters from the measurement set.
"""
function collect_metadata(ms::Table, beam)
    antennas  = read_antennas(ms)
    baselines = read_baselines(ms)
    channels  = read_channels(ms)
    phase_center = read_phase_center(ms)
    time = read_time(ms)
    Metadata(antennas, baselines, channels, phase_center, time, beam)
end

Nant(meta::Metadata)  = length(meta.antennas)
Nfreq(meta::Metadata) = length(meta.channels)
Nbase(meta::Metadata) = length(meta.baselines)

function position(meta)
    x = 0.0
    y = 0.0
    z = 0.0
    for i = 1:Nant(meta)
        pos = meta.antennas[i].position
        x += pos.x
        y += pos.y
        z += pos.z
    end
    x /= Nant(meta)
    y /= Nant(meta)
    z /= Nant(meta)
    Position(pos"ITRF", x, y, z)
end

function reference_frame(meta)
    frame = ReferenceFrame()
    set!(frame, meta.time)
    set!(frame, position(meta))
    frame
end

type Visibilities
    data  :: Matrix{JonesMatrix}
    flags :: Matrix{Bool}
end

Nfreq(vis::Visibilities) = size(vis.data, 2)
Nbase(vis::Visibilities) = size(vis.data, 1)

"""
    read_data_column(ms::Table, column)

Read the visibilities from the measurement set.
"""
function read_data_column(ms::Table, column)
    raw_data   = ms[column]
    data_flags = ms["FLAG"]
    row_flags  = ms["FLAG_ROW"]
    data  = organize_data(raw_data)
    flags = resolve_flags(data_flags, row_flags)
    Visibilities(data, flags)
end

function write_data_column(ms::Table, column, data::Visibilities)
    reorganized_data = zeros(Complex64, 4, Nfreq(data), Nbase(data))
    for α = 1:Nbase(data), β = 1:Nfreq(data)
        reorganized_data[1,β,α] = data.data[α,β].xx
        reorganized_data[2,β,α] = data.data[α,β].xy
        reorganized_data[3,β,α] = data.data[α,β].yx
        reorganized_data[4,β,α] = data.data[α,β].yy
    end
    ms[column] = reorganized_data
end

"""
    read_antennas(ms::Table)

Read the antenna positions from the `ANTENNA` subtable.
"""
function read_antennas(ms::Table)
    antenna_table = ms[kw"ANTENNA"] |> Table
    xyz = antenna_table["POSITION"]
    antennas = Antenna[]
    for i = 1:size(xyz, 2)
        x = xyz[1,i]
        y = xyz[2,i]
        z = xyz[3,i]
        position = Position(pos"ITRF", x, y, z)
        push!(antennas, Antenna(position))
    end
    unlock(antenna_table)
    antennas
end

function read_baselines(ms::Table)
    ant1 = ms["ANTENNA1"]
    ant2 = ms["ANTENNA2"]
    [Baseline(ant1[α]+1, ant2[α]+1) for α = 1:length(ant1)]
end

function read_channels(ms::Table)
    spw_table = ms[kw"SPECTRAL_WINDOW"] |> Table
    channels  = spw_table["CHAN_FREQ", 1]
    unlock(spw_table)
    channels
end

function read_phase_center(ms::Table)
    field_table = ms[kw"FIELD"] |> Table
    dir = field_table["PHASE_DIR"]
    unlock(field_table)
    Direction(dir"J2000", dir[1]*radians, dir[2]*radians)
end

function read_time(ms::Table)
    time = ms["TIME", 1]
    Epoch(epoch"UTC", time*seconds)
end

function organize_data(raw_data)
    data = zeros(JonesMatrix, size(raw_data,3), size(raw_data,2))
    for α = 1:size(raw_data,3), β = 1:size(raw_data,2)
        data[α,β] = JonesMatrix(raw_data[1,β,α], raw_data[2,β,α],
                                raw_data[3,β,α], raw_data[4,β,α])
    end
    data
end

function resolve_flags(data_flags, row_flags)
    flags = zeros(Bool, size(data_flags,3), size(data_flags,2))
    for β = 1:size(data_flags,2), α = 1:size(data_flags,3)
        if data_flags[1,β,α] || data_flags[2,β,α] || data_flags[3,β,α] || data_flags[4,β,α]
            flags[α,β] = true
        end
    end
    for α = 1:length(row_flags)
        if row_flags[α]
            flags[α,:] = true
        end
    end
    flags
end

function get_data(ms::Table)
    read_data_column(ms, "DATA")
end

function set_data!(ms::Table, data)
    write_data_column(ms, "DATA", data)
end

function get_model_data(ms::Table)
    read_data_column(ms, "MODEL_DATA")
end

function set_model_data!(ms::Table, data, force_imaging_columns = false)
    if force_imaging_columns || Tables.exists(ms,"MODEL_DATA")
        write_data_column(ms, "MODEL_DATA", data)
    end
end

"""
    get_corrected_data(ms::MeasurementSet)

Get the CORRECTED_DATA column if it exists. Otherwise settle for the DATA column.
"""
function get_corrected_data(ms::Table)
    if Tables.exists(ms, "CORRECTED_DATA")
        return read_data_column(ms, "CORRECTED_DATA")
    else
        return read_data_column(ms, "DATA")
    end
end

function set_corrected_data!(ms::Table, data, force_imaging_columns = false)
    if force_imaging_columns || Tables.exists(ms, "CORRECTED_DATA")
        write_data_column(ms, "CORRECTED_DATA", data)
    else
        write_data_column(ms, "DATA", data)
    end
end

#function set_flags!(ms::MeasurementSet,flags)
#    ms.table["FLAG"] = flags
#end

function flag_short_baselines!(data, meta, minuvw)
    for β = 1:Nfreq(meta)
        ν = meta.channels[β]
        λ = c / ν
        for α = 1:Nbase(meta)
            antenna1 = meta.antennas[meta.baselines[α].antenna1]
            antenna2 = meta.antennas[meta.baselines[α].antenna2]
            u = antenna1.position.x - antenna2.position.x
            v = antenna1.position.y - antenna2.position.y
            w = antenna1.position.z - antenna2.position.z
            b = sqrt(u^2 + v^2 + w^2)
            if b < minuvw * λ
                data.flags[α,β] = true
            end
        end
    end
    data
end

