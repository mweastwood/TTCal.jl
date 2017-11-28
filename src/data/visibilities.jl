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

abstract type Visibilities{T} end

struct SinglePolarizationVisibilities{T} <: Visibilities{T}
    data  :: Vector{ComplexVector{T}}
    flags :: Vector{Vector{Bool}}
end

function SinglePolarizationVisibilities(Nant)
    data  = [ComplexVector(Nant) for antenna = 1:Nant]
    flags = [   fill(true, Nant) for antenna = 1:Nant]
    SinglePolarizationVisibilities(data, flags)
end

#struct DualPolarizationVisibilities{T}
#end

#struct FullPolarizationVisibilities{T}
#end

function Base.getindex(visibilities::Visibilities, antenna1, antenna2)
    value = visibilities.data[antenna2][antenna1]
    flag  = visibilities.flags[antenna2][antenna1]
    value, flag
end

function Base.setindex!(visibilities::Visibilities, value, antenna1, antenna2)
    visibilities.data[antenna2][antenna1]  = value
    visibilities.data[antenna1][antenna2]  = value'
    visibilities.flags[antenna2][antenna1] = false
    visibilities.flags[antenna1][antenna2] = false
    value
end

function flag!(visibilities::Visibilities, antenna1, antenna2)
    visibilities.data[antenna2][antenna1]  = 0
    visibilities.flags[antenna2][antenna1] = true
    nothing
end

mutable struct Dataset{T <: Visibilities}
    metadata :: Metadata
    data     :: Matrix{T} # (number of frequency channels) × (number of time integrations)
end

function Dataset(T::Type{<:Visibilities}, metadata)
    data = [T(Nant(metadata)) for freq = 1:Nfreq(metadata), time = 1:Ntime(metadata)]
    Dataset(metadata, data)
end

function Base.getindex(dataset::Dataset, channel, time)
    dataset.data[channel, time]
end

Nfreq(dataset::Dataset) = Nfreq(dataset.metadata)
Ntime(dataset::Dataset) = Ntime(dataset.metadata)
Nant(dataset::Dataset)  =  Nant(dataset.metadata)
Nbase(dataset::Dataset) = Nbase(dataset.metadata)

function merge!(lhs::Dataset{T}, rhs::Dataset{T}; axis=:frequency) where {T}
    merge!(lhs.metadata, rhs.metadata, axis=axis)
    Nfreq,  Ntime  = size(lhs.data)
    Nfreq′, Ntime′ = size(rhs.data)
    if axis == :frequency
        data = Array{T}(Nfreq+Nfreq′, Ntime)
        data[    1:Nfreq, :] = lhs.data
        data[Nfreq+1:end, :] = rhs.data
        lhs.data = data
    elseif axis == :time
        data = Array{T}(Nfreq, Ntime+Ntime′)
        data[:,     1:Ntime] = lhs.data
        data[:, Ntime+1:end] = rhs.data
        lhs.data = data
    else
        err("unknown merge axis $axis")
    end
    lhs
end

function slice!(dataset::Dataset{T}, indices; axis=:frequency) where {T}
    slice!(dataset.metadata, indices, axis=axis)
    if axis == :frequency
        dataset.data = dataset.data[indices, :]
    elseif axis == :time
        dataset.data = dataset.data[:, indices]
    else
        err("unknown slice axis $axis")
    end
end






"Read dataset from the given measurement set."
function Dataset(ms::Table; column="DATA", polarization=:full)
    metadata = Metadata(ms::Table)
    if polarization == :aa
        return read_single_polarization_visibilities(ms, metadata, column, 1)
    elseif polarization == :bb
        return read_single_polarization_visibilities(ms, metadata, column, 4)
    else
        err("unknown polarization $polarization")
    end
end

function read_single_polarization_visibilities(ms, metadata, column, polarization)
    dataset = Dataset(SinglePolarizationVisibilities, metadata)
    data  = ms[column]
    flags = read_flags(ms)
    antenna1, antenna2 = read_antenna1_antenna2(ms)
    for frequency = 1:Nfreq(metadata)
        visibilities = dataset[frequency, 1]
        for baseline = 1:Nbase(metadata)
            flags[baseline, frequency] && continue
            ant1 = antenna1[baseline]
            ant2 = antenna2[baseline]
            visibilities[ant1, ant2] = data[polarization, frequency, baseline]
        end
    end
    dataset
end

"Read flags from the given measurement set."
function read_flags(ms::Table)
    data_flags = ms["FLAG"]
    row_flags  = ms["FLAG_ROW"]
    resolve_flags(data_flags, row_flags)
end

function resolve_flags(data_flags, row_flags)
    # Merge the row flags with the rest of the flags.
    # Note that here we flag all polarizations if one polarization is flagged.
    # (we may want to change this behavior in the future)
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

"Read the ANTENNA1 and ANTENNA2 columns from the given measurement set."
function read_antenna1_antenna2(ms::Table)
    # Note that we add 1 so that these antenna numbers start at 1
    antenna1 = ms["ANTENNA1"] .+ 1
    antenna2 = ms["ANTENNA2"] .+ 1
    antenna1, antenna2
end








#"Write visibilities to the measurement set."
#function write(ms::Table, column, data::Visibilities; apply_flags::Bool=true)
#    reorganized_data = zeros(Complex64, 4, Nfreq(data), Nbase(data))
#    for α = 1:Nbase(data), β = 1:Nfreq(data)
#        reorganized_data[1,β,α] = data.data[α,β].xx
#        reorganized_data[2,β,α] = data.data[α,β].xy
#        reorganized_data[3,β,α] = data.data[α,β].yx
#        reorganized_data[4,β,α] = data.data[α,β].yy
#    end
#    ms[column] = reorganized_data
#    if apply_flags
#        flags = zeros(Bool, 4, Nfreq(data), Nbase(data))
#        for α = 1:Nbase(data), β = 1:Nfreq(data)
#            if data.flags[α,β]
#                flags[:,β,α] = true
#            end
#        end
#        ms["FLAG"] = flags
#    end
#end
#
#function organize_data(raw_data)
#    data = zeros(JonesMatrix, size(raw_data,3), size(raw_data,2))
#    for α = 1:size(raw_data,3), β = 1:size(raw_data,2)
#        data[α,β] = JonesMatrix(raw_data[1,β,α], raw_data[2,β,α],
#                                raw_data[3,β,α], raw_data[4,β,α])
#    end
#    data
#end
#

#function flag_short_baselines!(data, meta, minuvw)
#    for β = 1:Nfreq(meta)
#        ν = meta.channels[β]
#        λ = c / ν
#        for α = 1:Nbase(meta)
#            antenna1 = meta.antennas[meta.baselines[α].antenna1]
#            antenna2 = meta.antennas[meta.baselines[α].antenna2]
#            u = antenna1.position.x - antenna2.position.x
#            v = antenna1.position.y - antenna2.position.y
#            w = antenna1.position.z - antenna2.position.z
#            b = sqrt(u^2 + v^2 + w^2)
#            if b < minuvw * λ
#                data.flags[α,β] = true
#            end
#        end
#    end
#    data
#end

