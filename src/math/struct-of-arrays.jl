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

"Implement the \"struct of arrays\" optimization to facilitiate SIMD vectorization."
struct ComplexVector{T} <: AbstractVector{Complex{T}}
    real :: Vector{T}
    imag :: Vector{T}
end

ComplexVector(N) = ComplexVector(Complex128, N)
ComplexVector(T, N) = ComplexVector(zeros(real(T), N), zeros(real(T), N))
Base.size(v::ComplexVector) = size(v.real)
Base.length(v::ComplexVector) = length(v.real)
Base.eltype(v::ComplexVector{T}) where {T} = complex(T)
Base.@propagate_inbounds function Base.getindex(v::ComplexVector, i)
    complex(v.real[i], v.imag[i])
end
Base.@propagate_inbounds function Base.setindex!(v::ComplexVector, value, i)
    v.real[i] = real(value)
    v.imag[i] = imag(value)
    value
end

