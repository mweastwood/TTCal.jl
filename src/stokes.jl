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

doc"""
    StokesVector

This type represents a Stokes vector.

``` math
    \begin{pmatrix}
        I \\
        Q \\
        U \\
        V \\
    \end{pmatrix}
```
"""
immutable StokesVector
    I::Float64
    Q::Float64
    U::Float64
    V::Float64
end

StokesVector() = zero(StokesVector)
Base.zero(::Type{StokesVector}) = StokesVector(0, 0, 0, 0)
Base.one(::Type{StokesVector}) = StokesVector(1, 0, 0, 0)
Base.rand(::Type{StokesVector}) = StokesVector(rand(), rand(), rand(), rand())

function StokesVector(v::Vector)
    length(v) == 4 || throw(DimensionMismatch("A Stokes vector must have 4 elements."))
    StokesVector(v[1],v[2],v[3],v[4])
end

Base.convert(::Type{Vector{Float64}}, v::StokesVector) = [v.I, v.Q, v.U, v.V]
Base.convert(::Type{Vector}, v::StokesVector) = Base.convert(Vector{Float64}, v)

function Base.show(io::IO, v::StokesVector)
    @printf(io, "(%.3f, %.3f, %.3f, %.3f)", v.I, v.Q, v.U, v.V)
end

function *(a::Number, v::StokesVector)
    StokesVector(a*v.I, a*v.Q, a*v.U, a*v.V)
end
*(v::StokesVector, a::Number) = a*v
/(v::StokesVector, a::Number) = (1/a)*v

for op in (:+, :-)
    @eval function $op(v1::StokesVector, v2::StokesVector)
        StokesVector($op(v1.I,v2.I), $op(v1.Q,v2.Q), $op(v1.U,v2.U), $op(v1.V,v2.V))
    end
end

function Base.norm(v::StokesVector)
    sqrt(abs2(v.I)+abs2(v.Q)+abs2(v.U)+abs2(v.V))
end

"""
    MuellerMatrix

This type represents a Mueller matrix.

    MuellerMatrix(J::JonesMatrix)

Create a Mueller matrix from the given Jones matrix.
"""
immutable MuellerMatrix
    mat::Matrix{Float64}
end

function MuellerMatrix(J::JonesMatrix)
    MuellerMatrix(to_stokes*kron(J,conj(J))*to_linear |> real)
end

Base.convert(::Type{Matrix{Float64}}, M::MuellerMatrix) = M.mat
Base.convert(::Type{Matrix}, M::MuellerMatrix) = Base.convert(Matrix{Float64}, M)

for op in (:+, :-, :*, :/)
    @eval function $op(M1::MuellerMatrix, M2::MuellerMatrix)
        MuellerMatrix($op(M1.mat, M2.mat))
    end
end

function *(M::MuellerMatrix, v::StokesVector)
    Matrix(M)*Vector(v) |> StokesVector
end

Base.norm(M::MuellerMatrix) = vecnorm(M.mat)

# Note the factor of 0.5 appears to be a convention in radio astronomy (but not physics)
const to_stokes = [0.5+0.0im  0.0+0.0im  0.0+0.0im   0.5+0.0im
                   0.5+0.0im  0.0+0.0im  0.0+0.0im  -0.5+0.0im
                   0.0+0.0im  0.5+0.0im  0.5+0.0im   0.0+0.0im
                   0.0+0.0im  0.0+0.5im  0.0-0.5im   0.0+0.0im]
const to_linear = [1.0+0.0im   1.0+0.0im  0.0+0.0im  0.0+0.0im
                   0.0+0.0im   0.0+0.0im  1.0+0.0im  0.0-1.0im
                   0.0+0.0im   0.0+0.0im  1.0+0.0im  0.0+1.0im
                   1.0+0.0im  -1.0-0.0im  0.0+0.0im  0.0+0.0im]

function StokesVector(correlations::HermitianJonesMatrix)::StokesVector
    vec = [correlations.xx, correlations.xy, conj(correlations.xy), correlations.yy]
    to_stokes*vec |> real |> StokesVector
end

function HermitianJonesMatrix(stokes::StokesVector)::HermitianJonesMatrix
    correlations = to_linear*Vector(stokes)
    HermitianJonesMatrix(real(correlations[1]), correlations[2], real(correlations[4]))
end

