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

struct Visibilities{P <: Polarization, T}
    data  :: Vector{T}
    flags :: Vector{Vector{Bool}}
end

function Visibilities(pol::Type{<:Polarization}, Nant)
    T = vector_type(pol)
    data  = [         T(Nant) for antenna = 1:Nant]
    flags = [fill(true, Nant) for antenna = 1:Nant]
    Visibilities{pol, T}(data, flags)
end

Nant( visibilities::Visibilities) = length(visibilities.data)
Nbase(visibilities::Visibilities) = Nbase(Nant(visibilities))
polarization(::Visibilities{P, T}) where {P, T} = P
Base.eltype( ::Visibilities{P, T}) where {P, T} = eltype(T)

Base.@propagate_inbounds function Base.getindex(visibilities::Visibilities, antenna)
    visibilities.data[antenna]
end

Base.@propagate_inbounds function Base.getindex(visibilities::Visibilities, antenna1, antenna2)
    visibilities.data[antenna2][antenna1]
end

Base.@propagate_inbounds function Base.setindex!(visibilities::Visibilities, value, ant1, ant2)
    visibilities.data[ant2][ant1]  = value
    visibilities.data[ant1][ant2]  = value'
    visibilities.flags[ant2][ant1] = false
    visibilities.flags[ant1][ant2] = false
    value
end

function isflagged(visibilities::Visibilities, antenna1, antenna2)
    visibilities.flags[antenna2][antenna1]
end

function flag!(visibilities::Visibilities, antenna1, antenna2)
    visibilities.data[antenna2][antenna1]  = zero(eltype(visibilities))
    visibilities.data[antenna1][antenna2]  = zero(eltype(visibilities))
    visibilities.flags[antenna2][antenna1] = true
    visibilities.flags[antenna1][antenna2] = true
    nothing
end

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

function Base.getindex(dataset::Dataset, channel, time)
    dataset.data[channel, time]
end

function flag_antennas!(dataset::Dataset, antennas)
    for time = 1:Ntime(dataset), frequency = 1:Nfreq(dataset)
        visibilities = dataset[frequency, time]
        for antenna1 in antennas, antenna2 = 1:Nant(dataset)
            flag!(visibilities, antenna1, antenna2)
        end
    end
    antennas
end

function flag_baselines!(dataset::Dataset, baselines)
    for time = 1:Ntime(dataset), frequency = 1:Nfreq(dataset)
        visibilities = dataset[frequency, time]
        for (antenna1, antenna2) in baselines
            flag!(visibilities, antenna1, antenna2)
        end
    end
    baselines
end

function flag_frequencies!(dataset::Dataset, frequencies)
    for time = 1:Ntime(dataset), frequency in frequencies
        visibilities = dataset[frequency, time]
        for antenna1 = 1:Nant(dataset), antenna2 = antenna1:Nant(dataset)
            flag!(visibilities, antenna1, antenna2)
        end
    end
    frequencies
end

function match_flags!(to, from)
    for time = 1:Ntime(to), frequency = 1:Nfreq(to)
        to_visibilities   =   to[frequency, time]
        from_visibilities = from[frequency, time]
        for antenna1 = 1:Nant(to), antenna2 = antenna1:Nant(to)
            if isflagged(from_visibilities, antenna1, antenna2)
                flag!(to_visibilities, antenna1, antenna2)
            end
        end
    end
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
    slice!(metadata, index:index, axis=axis)
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

#function read_dataset_jld2(file)
#    #metadata = file["METADATA"]
#    #pol    = file["POLARIZATION"]
#    #vector = file["DATA"]
#    #flags  = file["FLAGS"]
#    #read_dataset_jld2(pol, metadata, vector, flags)
#end

#function read_dataset_jld2(::Type{TTCal.Dual}, metadata, vector, flags)
#    dataset = Dataset(metadata, polarization=TTCal.Dual)
#    array = reshape(vector, (4, Nbase(dataset), Nfreq(dataset), Ntime(dataset)))
#    for time = 1:Ntime(dataset), frequency = 1:Nfreq(dataset)
#        visibilities = dataset[frequency, time]
#        baseline = 1
#        for antenna1 = 1:Nant(dataset), antenna2 = antenna1:Nant(dataset)
#            if !flags[baseline, frequency, time]
#                J = DiagonalJonesMatrix(complex(array[1, baseline, frequency, time],
#                                                array[2, baseline, frequency, time]),
#                                        complex(array[3, baseline, frequency, time],
#                                                array[4, baseline, frequency, time]))
#                visibilities[antenna1, antenna2] = J
#            end
#            baseline += 1
#        end
#    end
#    dataset
#end

#function write_dataset_jld2!(file, dataset::Dataset, ::Type{TTCal.Dual})
#    T = Float64
#    array = zeros(   T, 4, Nbase(dataset), Nfreq(dataset), Ntime(dataset))
#    flags = zeros(Bool,    Nbase(dataset), Nfreq(dataset), Ntime(dataset))
#    for time = 1:Ntime(dataset), frequency = 1:Nfreq(dataset)
#        visibilities = dataset[frequency, time]
#        baseline = 1
#        for antenna1 = 1:Nant(dataset), antenna2 = antenna1:Nant(dataset)
#            J = visibilities[antenna1, antenna2]
#            array[1, baseline, frequency, time] = real(J.xx)
#            array[2, baseline, frequency, time] = imag(J.xx)
#            array[3, baseline, frequency, time] = real(J.yy)
#            array[4, baseline, frequency, time] = imag(J.yy)
#            flags[   baseline, frequency, time] = isflagged(visibilities, antenna1, antenna2)
#            baseline += 1
#        end
#    end
#    file["DATA"]  = vec(array)
#    file["FLAGS"] = vec(flags)
#    file["POLARIZATION"] = TTCal.Dual
#end







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

