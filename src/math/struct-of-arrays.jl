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
ComplexVector{T}(N) where {T} = ComplexVector(T, N)
Base.size(v::ComplexVector) = size(v.real)
Base.length(v::ComplexVector) = length(v.real)
Base.eltype(v::ComplexVector{T}) where {T} = complex(T)
Base.@propagate_inbounds function Base.getindex(v::ComplexVector, i::Integer)
    complex(v.real[i], v.imag[i])
end
Base.@propagate_inbounds function Base.setindex!(v::ComplexVector, value, i)
    v.real[i] = real(value)
    v.imag[i] = imag(value)
    value
end

"Implement the \"struct of arrays\" optimization to facilitiate SIMD vectorization."
struct DiagonalJonesVector <: AbstractVector{DiagonalJonesMatrix}
    xxreal :: Vector{Float64}
    xximag :: Vector{Float64}
    yyreal :: Vector{Float64}
    yyimag :: Vector{Float64}
end

DiagonalJonesVector(N) = DiagonalJonesVector(zeros(Float64, N), zeros(Float64, N),
                                             zeros(Float64, N), zeros(Float64, N))
Base.size(v::DiagonalJonesVector) = size(v.xxreal)
Base.length(v::DiagonalJonesVector) = length(v.xxreal)
Base.eltype(v::DiagonalJonesVector) = DiagonalJonesMatrix
Base.@propagate_inbounds function Base.getindex(v::DiagonalJonesVector, i::Integer)
    DiagonalJonesMatrix(complex(v.xxreal[i], v.xximag[i]),
                        complex(v.yyreal[i], v.yyimag[i]))
end
Base.@propagate_inbounds function Base.setindex!(v::DiagonalJonesVector, value, i)
    v.xxreal[i] = real(value.xx)
    v.xximag[i] = imag(value.xx)
    v.yyreal[i] = real(value.yy)
    v.yyimag[i] = imag(value.yy)
    value
end

"Implement the \"struct of arrays\" optimization to facilitiate SIMD vectorization."
struct JonesVector <: AbstractVector{JonesMatrix}
    xxreal :: Vector{Float64}
    xximag :: Vector{Float64}
    xyreal :: Vector{Float64}
    xyimag :: Vector{Float64}
    yxreal :: Vector{Float64}
    yximag :: Vector{Float64}
    yyreal :: Vector{Float64}
    yyimag :: Vector{Float64}
end

JonesVector(N) = JonesVector(zeros(Float64, N), zeros(Float64, N),
                             zeros(Float64, N), zeros(Float64, N),
                             zeros(Float64, N), zeros(Float64, N),
                             zeros(Float64, N), zeros(Float64, N))
Base.size(v::JonesVector) = size(v.xxreal)
Base.length(v::JonesVector) = length(v.xxreal)
Base.eltype(v::JonesVector) = JonesMatrix
Base.@propagate_inbounds function Base.getindex(v::JonesVector, i::Integer)
    JonesMatrix(complex(v.xxreal[i], v.xximag[i]),
                complex(v.xyreal[i], v.xyimag[i]),
                complex(v.yxreal[i], v.yximag[i]),
                complex(v.yyreal[i], v.yyimag[i]))
end
Base.@propagate_inbounds function Base.setindex!(v::JonesVector, value, i)
    v.xxreal[i] = real(value.xx)
    v.xximag[i] = imag(value.xx)
    v.xyreal[i] = real(value.xy)
    v.xyimag[i] = imag(value.xy)
    v.yxreal[i] = real(value.yx)
    v.yximag[i] = imag(value.yx)
    v.yyreal[i] = real(value.yy)
    v.yyimag[i] = imag(value.yy)
    value
end

vector_type(::Type{XX})   = ComplexVector{Float64}
vector_type(::Type{XY})   = ComplexVector{Float64}
vector_type(::Type{YX})   = ComplexVector{Float64}
vector_type(::Type{YY})   = ComplexVector{Float64}
vector_type(::Type{Dual}) = DiagonalJonesVector
vector_type(::Type{Full}) = JonesVector

