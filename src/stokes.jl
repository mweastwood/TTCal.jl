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

doc"""
    immutable JonesMatrix

This type represents a 2x2 complex Jones matrix.

\\[
\begin{pmatrix}
    j_{xx} & j_{xy} \\\\
    j_{yx} & j_{yy} \\\\
\end{pmatrix}
\\]
"""
immutable JonesMatrix
    xx::Complex128
    xy::Complex128
    yx::Complex128
    yy::Complex128
end

doc"""
    immutable DiagonalJonesMatrix

This type represents a Jones matrix that is diagonal.

These matrices are used to represent the complex gains of each antenna
without accounting for the off-diagonal polarization leakage terms.
"""
immutable DiagonalJonesMatrix
    xx::Complex128
    yy::Complex128
end

doc"""
    immutable HermitianJonesMatrix

This type represents a Jones matrix that is Hermitian.

These matrices are useful for representing the $xx$, $xy$, $yx$, and $yy$ flux,
because $xx$ and $yy$ are constrained to be real while $xy$ and $yx$ are complex
conjugates.
"""
immutable HermitianJonesMatrix
    xx::Float64
    xy::Complex128
    yy::Float64
end

zero(::Type{JonesMatrix}) = JonesMatrix(0,0,0,0)
one(::Type{JonesMatrix}) = JonesMatrix(1,0,0,1) # the identity matrix
rand(::Type{JonesMatrix}) = JonesMatrix(rand(Complex128),rand(Complex128),rand(Complex128),rand(Complex128))

zero(::Type{DiagonalJonesMatrix}) = DiagonalJonesMatrix(0,0)
one(::Type{DiagonalJonesMatrix}) = DiagonalJonesMatrix(1,1) # the identity matrix
rand(::Type{DiagonalJonesMatrix}) = DiagonalJonesMatrix(rand(Float64),rand(Float64))

zero(::Type{HermitianJonesMatrix}) = HermitianJonesMatrix(0,0,0)
one(::Type{HermitianJonesMatrix}) = HermitianJonesMatrix(1,0,1) # the identity matrix
rand(::Type{HermitianJonesMatrix}) = HermitianJonesMatrix(rand(Float64),rand(Complex128),rand(Float64))

function JonesMatrix(mat::Matrix)
    size(mat) == (2,2) || throw(DimensionMismatch("A Jones matrix must be 2x2."))
    JonesMatrix(mat[1,1],mat[1,2],mat[2,1],mat[2,2])
end

function Base.convert(::Type{Matrix{Complex128}},J::JonesMatrix)
    [J.xx J.xy
     J.yx J.yy]
end
Base.convert(::Type{Matrix},J::JonesMatrix) = Base.convert(Matrix{Complex128},J)

function Base.convert(::Type{Matrix{Complex128}},J::HermitianJonesMatrix)
    [J.xx  J.xy
     J.xy' J.yy]
end
Base.convert(::Type{Matrix},J::HermitianJonesMatrix) = Base.convert(Matrix{Complex128},J)

for op in (:+,:-)
    @eval @inline function $op(J1::JonesMatrix,J2::JonesMatrix)
        JonesMatrix($op(J1.xx,J2.xx),
                    $op(J1.xy,J2.xy),
                    $op(J1.yx,J2.yx),
                    $op(J1.yy,J2.yy))
    end
    @eval @inline function $op(J1::DiagonalJonesMatrix,J2::DiagonalJonesMatrix)
        DiagonalJonesMatrix($op(J1.xx,J2.xx),
                            $op(J1.yy,J2.yy))
    end
    @eval @inline function $op(J1::HermitianJonesMatrix,J2::HermitianJonesMatrix)
        HermitianJonesMatrix($op(J1.xx,J2.xx),
                             $op(J1.xy,J2.xy),
                             $op(J1.yy,J2.yy))
    end
end

@inline *(a::Number,J::JonesMatrix) = JonesMatrix(a*J.xx,a*J.xy,a*J.yx,a*J.yy)
@inline *(J::JonesMatrix,a::Number) = *(a,J)

@inline *(a::Number,J::DiagonalJonesMatrix) = DiagonalJonesMatrix(a*J.xx,a*J.yy)
@inline *(J::DiagonalJonesMatrix,a::Number) = *(a,J)

@inline /(J::HermitianJonesMatrix,a::Number) = HermitianJonesMatrix(J.xx/a,J.xy/a,J.yy/a)

@inline function *(J1::JonesMatrix,J2::JonesMatrix)
    JonesMatrix(J1.xx*J2.xx + J1.xy*J2.yx,
                J1.xx*J2.xy + J1.xy*J2.yy,
                J1.yx*J2.xx + J1.yy*J2.yx,
                J1.yx*J2.xy + J1.yy*J2.yy)
end

@inline function *(J1::JonesMatrix,J2::DiagonalJonesMatrix)
    JonesMatrix(J1.xx*J2.xx,J1.xy*J2.yy,
                J1.yx*J2.xx,J1.yy*J2.yy)
end

@inline function *(J1::DiagonalJonesMatrix,J2::JonesMatrix)
    JonesMatrix(J1.xx*J2.xx,J1.xx*J2.xy,
                J1.yy*J2.yx,J1.yy*J2.yy)
end

@inline *(J1::DiagonalJonesMatrix,J2::DiagonalJonesMatrix) = DiagonalJonesMatrix(J1.xx*J2.xx,J1.yy*J2.yy)

@inline \(J1::JonesMatrix,J2::JonesMatrix) = inv(J1)*J2
@inline \(J1::DiagonalJonesMatrix,J2::DiagonalJonesMatrix) = DiagonalJonesMatrix(J2.xx/J1.xx,J2.yy/J1.yy)

@inline conj(J::JonesMatrix) = JonesMatrix(conj(J.xx),conj(J.xy),conj(J.yx),conj(J.yy))
@inline conj(J::DiagonalJonesMatrix) = DiagonalJonesMatrix(conj(J.xx),conj(J.yy))

@inline ctranspose(J::JonesMatrix) = JonesMatrix(conj(J.xx),conj(J.yx),conj(J.xy),conj(J.yy))
@inline ctranspose(J::DiagonalJonesMatrix) = conj(J)

@inline det(J::JonesMatrix) = J.xx*J.yy - J.xy*J.yx
@inline det(J::DiagonalJonesMatrix) = J.xx*J.yy

function inv(J::JonesMatrix)
    d = det(J)
    JonesMatrix(J.yy/d,-J.xy/d,-J.yx/d,J.xx/d)
end

function inv(J::DiagonalJonesMatrix)
    DiagonalJonesMatrix(1/J.xx,1/J.yy)
end

# use the Frobenius norm
norm(J::JonesMatrix) = sqrt(abs2(J.xx)+abs2(J.xy)+abs2(J.yx)+abs2(J.yy))
norm(J::DiagonalJonesMatrix) = sqrt(abs2(J.xx)+abs2(J.yy))

function kron(J1::JonesMatrix,J2::JonesMatrix)
    [J1.xx*J2.xx J1.xx*J2.xy J1.xy*J2.xx J1.xy*J2.xy;
     J1.xx*J2.yx J1.xx*J2.yy J1.xy*J2.yx J1.xy*J2.yy;
     J1.yx*J2.xx J1.yx*J2.xy J1.yy*J2.xx J1.yy*J2.xy;
     J1.yx*J2.yx J1.yx*J2.yy J1.yy*J2.yx J1.yy*J2.yy]
end

doc"""
    congruence_transform(J::JonesMatrix, K::HermitianJonesMatrix) -> J*K*J'

Compute the congruence transformation of $K$ with respect to $J$:

\\[
    K \rightarrow JKJ^*
\\]
"""
function congruence_transform(J::JonesMatrix,K::HermitianJonesMatrix)
    HermitianJonesMatrix(abs2(J.xx)*K.xx + 2real(J.xx*J.xy'*K.xy) + abs2(J.xy)*K.yy,
                         J.xx*J.yx'*K.xx + J.xx*J.yy'*K.xy
                            + J.xy*J.yx'*K.xy' + J.xy*J.yy'*K.yy,
                         abs2(J.yx)*K.xx + 2real(J.yx*J.yy'*K.xy) + abs2(J.yy)*K.yy)
end

"""
    immutable MuellerMatrix

This type represents a Mueller matrix.
"""
immutable MuellerMatrix
    mat::Matrix{Float64}
end

"""
    MuellerMatrix(J::JonesMatrix)

Create a Mueller matrix from the given Jones matrix.
"""
function MuellerMatrix(J::JonesMatrix)
    MuellerMatrix(to_stokes*kron(J,conj(J))*to_linear |> real)
end

# Note the factor of 0.5 appears to be a convention in radio astronomy (but not physics)
const to_stokes = [0.5+0.0im  0.0+0.0im  0.0+0.0im   0.5+0.0im
                   0.5+0.0im  0.0+0.0im  0.0+0.0im  -0.5+0.0im
                   0.0+0.0im  0.5+0.0im  0.5+0.0im   0.0+0.0im
                   0.0+0.0im  0.0+0.5im  0.0-0.5im   0.0+0.0im]
const to_linear = [1.0+0.0im   1.0+0.0im  0.0+0.0im  0.0+0.0im
                   0.0+0.0im   0.0+0.0im  1.0+0.0im  0.0-1.0im
                   0.0+0.0im   0.0+0.0im  1.0+0.0im  0.0+1.0im
                   1.0+0.0im  -1.0-0.0im  0.0+0.0im  0.0+0.0im]

Base.convert(::Type{Matrix{Float64}},M::MuellerMatrix) = M.mat
Base.convert(::Type{Matrix},M::MuellerMatrix) = Base.convert(Matrix{Float64},M)

for op in (:+,:-,:*,:/)
    @eval function $op(M1::MuellerMatrix,M2::MuellerMatrix)
        MuellerMatrix($op(M1.mat,M2.mat))
    end
end

norm(M::MuellerMatrix) = vecnorm(M.mat)

doc"""
    immutable StokesVector

This type represents a Stokes vector.

\\[
\begin{pmatrix}
    I \\\\
    Q \\\\
    U \\\\
    V \\\\
\end{pmatrix}
\\]
"""
immutable StokesVector
    I::Float64
    Q::Float64
    U::Float64
    V::Float64
end

StokesVector() = zero(StokesVector)
zero(::Type{StokesVector}) = StokesVector(0,0,0,0)
rand(::Type{StokesVector}) = StokesVector(rand(),rand(),rand(),rand())

function StokesVector(v::Vector)
    length(v) == 4 || throw(DimensionMismatch("A Stokes vector must have 4 elements."))
    StokesVector(v[1],v[2],v[3],v[4])
end

Base.convert(::Type{Vector{Float64}},v::StokesVector) = [v.I,v.Q,v.U,v.V]
Base.convert(::Type{Vector},v::StokesVector) = Base.convert(Vector{Float64},v)

for op in (:+,:-)
    @eval function $op(v1::StokesVector,v2::StokesVector)
        StokesVector($op(v1.I,v2.I),
                     $op(v1.Q,v2.Q),
                     $op(v1.U,v2.U),
                     $op(v1.V,v2.V))
    end
end

function *(M::MuellerMatrix,v::StokesVector)
    Matrix(M)*Vector(v) |> StokesVector
end

function *(a::Number,v::StokesVector)
    StokesVector(a*v.I,a*v.Q,a*v.U,a*v.V)
end
*(v::StokesVector,a::Number) = *(a,v)

function norm(v::StokesVector)
    # the Frobenius norm
    sqrt(abs2(v.I)+abs2(v.Q)+abs2(v.U)+abs2(v.V))
end

doc"""
    stokes(correlations::HermitianJonesMatrix) -> StokesVector

Take the set of correlations $(xx,xy,yx,yy)$ and
convert it to the set of Stokes parameters $(I,Q,U,V)$.
"""
function stokes(correlations::HermitianJonesMatrix)
    vec = [correlations.xx, correlations.xy, conj(correlations.xy), correlations.yy]
    to_stokes*vec |> real |> StokesVector
end

doc"""
    linear(stokes::StokesVector) -> HermitianJonesMatrix

Take the set of Stokes parameters $(I,Q,U,V)$ and
convert it to the set of correlations $(xx,xy,yx,yy)$.
"""
function linear(stokes::StokesVector)
    correlations = to_linear*Vector(stokes)
    HermitianJonesMatrix(real(correlations[1]),
                              correlations[2],
                         real(correlations[4]))
end

