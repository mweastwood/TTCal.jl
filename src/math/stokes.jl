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

doc"""
    StokesVector

This type represents a Stokes vector.

```math
    \begin{pmatrix}
        I \\
        Q \\
        U \\
        V \\
    \end{pmatrix}
```
"""
immutable StokesVector <: AbstractVector{Float64}
    I::Float64
    Q::Float64
    U::Float64
    V::Float64
end

# Abstract array interface

Base.size(::StokesVector) = (4,)
Base.length(::StokesVector) = 4

function Base.getindex(v::StokesVector, idx::Int)
    @boundscheck checkbounds(v, idx)
    if idx == 1
        v.I
    elseif idx == 2
        v.Q
    elseif idx == 3
        v.U
    else
        v.V
    end
end

StokesVector() = zero(StokesVector)
Base.zero(::Type{StokesVector}) = StokesVector(0, 0, 0, 0)
Base.one( ::Type{StokesVector}) = StokesVector(1, 0, 0, 0)
Base.rand(::Type{StokesVector}) = StokesVector(rand(), rand(), rand(), rand())

for op in (:+, :-)
    @eval function Base.$op(v1::StokesVector, v2::StokesVector)
        StokesVector($op(v1.I, v2.I), $op(v1.Q, v2.Q),
                     $op(v1.U, v2.U), $op(v1.V, v2.V))
    end
end

for op in (:*, :/)
    @eval function Base.$op(v::StokesVector, a::Number)
        StokesVector($op(v.I, a), $op(v.Q, a), $op(v.U, a), $op(v.V, a))
    end
end
Base.:*(a::Number, v::StokesVector) = v*a

#function Base.norm(v::StokesVector)
#    sqrt(abs2(v.I)+abs2(v.Q)+abs2(v.U)+abs2(v.V))
#end

#"""
#    MuellerMatrix
#
#This type represents a Mueller matrix.
#
#    MuellerMatrix(J::JonesMatrix)
#
#Create a Mueller matrix from the given Jones matrix.
#"""
#immutable MuellerMatrix
#    mat::Matrix{Float64}
#end
#
#function MuellerMatrix(J::JonesMatrix)
#    MuellerMatrix(to_stokes*kron(J,conj(J))*to_linear |> real)
#end
#
#Base.convert(::Type{Matrix{Float64}}, M::MuellerMatrix) = M.mat
#Base.convert(::Type{Matrix}, M::MuellerMatrix) = Base.convert(Matrix{Float64}, M)
#
#for op in (:+, :-, :*, :/)
#    @eval function $op(M1::MuellerMatrix, M2::MuellerMatrix)
#        MuellerMatrix($op(M1.mat, M2.mat))
#    end
#end
#
#function *(M::MuellerMatrix, v::StokesVector)
#    Matrix(M)*Vector(v) |> StokesVector
#end
#
#Base.norm(M::MuellerMatrix) = vecnorm(M.mat)

# Note the factor of 0.5 appears to be a convention in radio astronomy (but not physics)
const to_stokes = @SMatrix [0.5+0.0im  0.0+0.0im  0.0+0.0im   0.5+0.0im
                            0.5+0.0im  0.0+0.0im  0.0+0.0im  -0.5+0.0im
                            0.0+0.0im  0.5+0.0im  0.5+0.0im   0.0+0.0im
                            0.0+0.0im  0.0+0.5im  0.0-0.5im   0.0+0.0im]
const to_linear = @SMatrix [1.0+0.0im   1.0+0.0im  0.0+0.0im  0.0+0.0im
                            0.0+0.0im   0.0+0.0im  1.0+0.0im  0.0-1.0im
                            0.0+0.0im   0.0+0.0im  1.0+0.0im  0.0+1.0im
                            1.0+0.0im  -1.0-0.0im  0.0+0.0im  0.0+0.0im]

function StokesVector(correlations::HermitianJonesMatrix)::StokesVector
    vec = @SVector [correlations.xx, correlations.xy, conj(correlations.xy), correlations.yy]
    vec = to_stokes * vec
    StokesVector(real(vec[1]), real(vec[2]), real(vec[3]), real(vec[4]))
end

function HermitianJonesMatrix(stokes::StokesVector)::HermitianJonesMatrix
    vec = to_linear*stokes
    HermitianJonesMatrix(real(vec[1]), vec[2], real(vec[4]))
end

