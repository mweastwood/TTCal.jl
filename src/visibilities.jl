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

type Visibilities
    data  :: Matrix{JonesMatrix}
    flags :: Matrix{Bool}
end

function Visibilities(Nbase::Int, Nfreq::Int)
    data = zeros(JonesMatrix, Nbase, Nfreq)
    flags = fill(false, Nbase, Nfreq)
    Visibilities(data, flags)
end

function Visibilities(meta::Metadata)
    Visibilities(Nbase(meta), Nfreq(meta))
end

Nbase(vis::Visibilities) = size(vis.data, 1)
Nfreq(vis::Visibilities) = size(vis.data, 2)

"Read visibilities from the measurement set."
function read(ms::Table, column)
    raw_data   = ms[column]
    data_flags = ms["FLAG"]
    row_flags  = ms["FLAG_ROW"]
    data  = organize_data(raw_data)
    flags = resolve_flags(data_flags, row_flags)
    Visibilities(data, flags)
end

"Write visibilities to the measurement set."
function write(ms::Table, column, data::Visibilities; apply_flags::Bool=true)
    reorganized_data = zeros(Complex64, 4, Nfreq(data), Nbase(data))
    for α = 1:Nbase(data), β = 1:Nfreq(data)
        reorganized_data[1,β,α] = data.data[α,β].xx
        reorganized_data[2,β,α] = data.data[α,β].xy
        reorganized_data[3,β,α] = data.data[α,β].yx
        reorganized_data[4,β,α] = data.data[α,β].yy
    end
    ms[column] = reorganized_data
    if apply_flags
        flags = zeros(Bool, 4, Nfreq(data), Nbase(data))
        for α = 1:Nbase(data), β = 1:Nfreq(data)
            if data.flags[α,β]
                flags[:,β,α] = true
            end
        end
        ms["FLAG"] = flags
    end
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

