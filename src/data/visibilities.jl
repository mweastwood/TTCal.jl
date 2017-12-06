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

