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

mutable struct Dataset{T <: Visibilities}
    metadata :: Metadata
    data     :: Matrix{T} # (number of frequency channels) × (number of time integrations)
end

function Dataset(metadata; polarization=Full)
    data = [Visibilities(polarization, Nant(metadata))
                for freq = 1:Nfreq(metadata), time = 1:Ntime(metadata)]
    Dataset(metadata, data)
end

function Base.rand!(dataset::Dataset)
    for time = 1:Ntime(dataset), frequency = 1:Nfreq(dataset)
        visibilities = dataset[frequency, time]
        for antenna1 = 1:Nant(dataset), antenna2 = antenna1:Nant(dataset)
            visibilities[antenna1, antenna2] = rand(eltype(dataset))
        end
    end
    dataset
end

function Base.getindex(dataset::Dataset, index)
    dataset.data[index]
end

function Base.getindex(dataset::Dataset, channel, time)
    dataset.data[channel, time]
end

Nfreq(dataset::Dataset) = Nfreq(dataset.metadata)
Ntime(dataset::Dataset) = Ntime(dataset.metadata)
Nant( dataset::Dataset) =  Nant(dataset.metadata)
Nbase(dataset::Dataset) = Nbase(dataset.metadata)
polarization(dataset::Dataset) = polarization(first(dataset.data))
Base.eltype( dataset::Dataset) = eltype(first(dataset.data))

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

function slice!(dataset::Dataset{T}, index::Integer; axis=:frequency) where{T}
    slice!(dataset, index:index, axis=axis)
end

function slice!(dataset::Dataset{T}, indices::AbstractVector; axis=:frequency) where {T}
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
function Dataset(ms::Table; column="DATA", polarization=Full)
    metadata = Metadata(ms)
    return read_dataset(ms, metadata, column, polarization)
end

function read_dataset(ms, metadata, column, polarization)
    dataset = Dataset(metadata, polarization=polarization)
    data  = ms[column]
    flags = read_flags(ms)
    antenna1, antenna2 = read_antenna1_antenna2(ms)
    pack_data!(dataset, data, flags, antenna1, antenna2, polarization)
    dataset
end

function pack_data!(dataset, data, flags, antenna1, antenna2, polarization)
    for frequency = 1:Nfreq(dataset)
        visibilities = dataset[frequency, 1]
        for baseline = 1:Nbase(dataset)
            flags[baseline, frequency] && continue
            ant1 = antenna1[baseline]
            ant2 = antenna2[baseline]
            element = read_dataset_element(data, frequency, baseline, polarization)
            visibilities[ant1, ant2] = element
        end
    end
    dataset
end

@inline polarization_index(::Type{XX}) = 1
@inline polarization_index(::Type{XY}) = 2
@inline polarization_index(::Type{YX}) = 3
@inline polarization_index(::Type{YY}) = 4

@inline function read_dataset_element(data, frequency, baseline, polarization::Type{<:Single})
    data[polarization_index(polarization), frequency, baseline]
end
@inline function read_dataset_element(data, frequency, baseline, polarization::Type{Dual})
    DiagonalJonesMatrix(read_dataset_element(data, frequency, baseline, XX),
                        read_dataset_element(data, frequency, baseline, YY))
end
@inline function read_dataset_element(data, frequency, baseline, polarization::Type{Full})
    JonesMatrix(read_dataset_element(data, frequency, baseline, XX),
                read_dataset_element(data, frequency, baseline, XY),
                read_dataset_element(data, frequency, baseline, YX),
                read_dataset_element(data, frequency, baseline, YY))
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

function write(ms::Table, dataset::Dataset; column="DATA")
    data  = zeros(Complex64, 4, Nfreq(dataset), Nbase(dataset))
    flags = zeros(     Bool, 4, Nfreq(dataset), Nbase(dataset))
    antenna1, antenna2 = read_antenna1_antenna2(ms)
    pol = polarization(dataset)
    for frequency = 1:Nfreq(dataset)
        visibilities = dataset[frequency, 1]
        for baseline = 1:Nbase(dataset)
            ant1 = antenna1[baseline]
            ant2 = antenna2[baseline]
            flags[:, frequency, baseline] = isflagged(visibilities, ant1, ant2)
            write_dataset_element!(data, frequency, baseline, visibilities[ant1, ant2], pol)
        end
    end
    if column == "DATA" || column == "CORRECTED_DATA"
        ms["FLAG"]     = flags
        ms["FLAG_ROW"] = zeros(Bool, Nbase(dataset))
    end
    ms[column] = data
end

function write_dataset_element!(data, frequency, baseline, element, polarization::Type{<:Single})
    data[polarization_index(polarization), frequency, baseline] = element
    element
end
function write_dataset_element!(data, frequency, baseline, element, polarization::Type{Dual})
    write_dataset_element!(data, frequency, baseline, element.xx, XX)
    write_dataset_element!(data, frequency, baseline, element.yy, YY)
    element
end
function write_dataset_element!(data, frequency, baseline, element, polarization::Type{Full})
    write_dataset_element!(data, frequency, baseline, element.xx, XX)
    write_dataset_element!(data, frequency, baseline, element.xy, XY)
    write_dataset_element!(data, frequency, baseline, element.yx, YX)
    write_dataset_element!(data, frequency, baseline, element.yy, YY)
    element
end

"Read dataset from the given JLD2 file."
function Dataset(path::String)
    local dataset
    jldopen(path, "r") do file
        dataset = file["DATASET"]
    end
    dataset
end

function write(path::String, dataset::Dataset)
    jldopen(path, "w", compress=true) do file
        file["DATASET"] = dataset
    end
end

