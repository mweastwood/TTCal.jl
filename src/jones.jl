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
    JonesMatrix

This type represents a 2x2 complex Jones matrix.

``` math
    \begin{pmatrix}
        j_{xx} & j_{xy} \\
        j_{yx} & j_{yy} \\
    \end{pmatrix}
```
"""
immutable JonesMatrix
    xx::Complex128
    xy::Complex128
    yx::Complex128
    yy::Complex128
end

doc"""
    DiagonalJonesMatrix

This type represents a Jones matrix that is diagonal.

``` math
    \begin{pmatrix}
        j_{xx} & 0 \\
        0 & j_{yy} \\
    \end{pmatrix}
```

These matrices are used to represent the complex gains of each antenna
without accounting for the off-diagonal polarization leakage terms.
"""
immutable DiagonalJonesMatrix
    xx::Complex128
    yy::Complex128
end

doc"""
    HermitianJonesMatrix

This type represents a Jones matrix that is Hermitian.

``` math
    \begin{pmatrix}
        j_{xx} & j_{xy} \\
        j_{xy}^* & j_{yy} \\
    \end{pmatrix}
```

These matrices are useful for representing the $xx$, $xy$, $yx$, and $yy$ flux,
because $xx$ and $yy$ are constrained to be real while $xy$ and $yx$ are complex
conjugates.
"""
immutable HermitianJonesMatrix
    xx::Float64
    xy::Complex128
    yy::Float64
end

typealias AnyJonesMatrix Union{JonesMatrix, DiagonalJonesMatrix, HermitianJonesMatrix}

Base.zero(::Type{JonesMatrix}) = JonesMatrix(0, 0, 0, 0)
Base.one(::Type{JonesMatrix}) = JonesMatrix(1, 0, 0, 1) # the identity matrix
Base.rand(::Type{JonesMatrix}) = JonesMatrix(rand(Complex128), rand(Complex128),
                                             rand(Complex128), rand(Complex128))

Base.zero(::Type{DiagonalJonesMatrix}) = DiagonalJonesMatrix(0, 0)
Base.one(::Type{DiagonalJonesMatrix}) = DiagonalJonesMatrix(1, 1) # the identity matrix
Base.rand(::Type{DiagonalJonesMatrix}) = DiagonalJonesMatrix(rand(Float64), rand(Float64))

Base.zero(::Type{HermitianJonesMatrix}) = HermitianJonesMatrix(0, 0, 0)
Base.one(::Type{HermitianJonesMatrix}) = HermitianJonesMatrix(1, 0, 1) # the identity matrix
Base.rand(::Type{HermitianJonesMatrix}) = HermitianJonesMatrix(rand(Float64), rand(Complex128),
                                                               rand(Float64))

function JonesMatrix(mat::Matrix)
    size(mat) == (2,2) || throw(DimensionMismatch("A Jones matrix must be 2x2."))
    JonesMatrix(mat[1,1], mat[1,2], mat[2,1], mat[2,2])
end

function Base.convert(::Type{Matrix{Complex128}}, J::JonesMatrix)
    [J.xx J.xy
     J.yx J.yy]
end
Base.convert(::Type{Matrix}, J::JonesMatrix) = Base.convert(Matrix{Complex128}, J)

function Base.convert(::Type{Matrix{Complex128}}, J::DiagonalJonesMatrix)
    [J.xx 0
     0    J.yy]
end
Base.convert(::Type{Matrix}, J::DiagonalJonesMatrix) = Base.convert(Matrix{Complex128}, J)

function Base.convert(::Type{Matrix{Complex128}}, J::HermitianJonesMatrix)
    [J.xx  J.xy
     J.xy' J.yy]
end
Base.convert(::Type{Matrix}, J::HermitianJonesMatrix) = Base.convert(Matrix{Complex128}, J)

for op in (:+, :-)
    @eval function $op(J1::JonesMatrix, J2::JonesMatrix)
        JonesMatrix($op(J1.xx, J2.xx),
                    $op(J1.xy, J2.xy),
                    $op(J1.yx, J2.yx),
                    $op(J1.yy, J2.yy))
    end
    @eval function $op(J1::DiagonalJonesMatrix, J2::DiagonalJonesMatrix)
        DiagonalJonesMatrix($op(J1.xx, J2.xx),
                            $op(J1.yy, J2.yy))
    end
    @eval function $op(J1::HermitianJonesMatrix, J2::HermitianJonesMatrix)
        HermitianJonesMatrix($op(J1.xx, J2.xx),
                             $op(J1.xy, J2.xy),
                             $op(J1.yy, J2.yy))
    end
end

*(a::Number, J::JonesMatrix) = JonesMatrix(a*J.xx, a*J.xy, a*J.yx, a*J.yy)
*(a::Number, J::DiagonalJonesMatrix) = DiagonalJonesMatrix(a*J.xx, a*J.yy)
*(a::Real, J::HermitianJonesMatrix) = HermitianJonesMatrix(a*J.xx, a*J.xy, a*J.yy)
*(a::Complex, J::HermitianJonesMatrix) = JonesMatrix(a*J.xx, a*J.xy, a*J.xy', a*J.yy)
*(J::AnyJonesMatrix, a::Number) = a*J

/(J::AnyJonesMatrix, a::Number) = (1/a)*J
./(J::AnyJonesMatrix, a::Number) = J/a

function *(J1::JonesMatrix, J2::JonesMatrix)
    JonesMatrix(J1.xx*J2.xx + J1.xy*J2.yx, J1.xx*J2.xy + J1.xy*J2.yy,
                J1.yx*J2.xx + J1.yy*J2.yx, J1.yx*J2.xy + J1.yy*J2.yy)
end

function *(J1::JonesMatrix, J2::DiagonalJonesMatrix)
    JonesMatrix(J1.xx*J2.xx,J1.xy*J2.yy, J1.yx*J2.xx,J1.yy*J2.yy)
end

function *(J1::DiagonalJonesMatrix, J2::JonesMatrix)
    JonesMatrix(J1.xx*J2.xx,J1.xx*J2.xy, J1.yy*J2.yx,J1.yy*J2.yy)
end

function *(J1::DiagonalJonesMatrix, J2::DiagonalJonesMatrix)
    DiagonalJonesMatrix(J1.xx*J2.xx, J1.yy*J2.yy)
end

function *(J1::HermitianJonesMatrix, J2::HermitianJonesMatrix)
    JonesMatrix(J1.xx *J2.xx + J1.xy*J2.xy', J1.xx *J2.xy + J1.xy*J2.yy,
                J1.xy'*J2.xx + J1.yy*J2.xy', J1.xy'*J2.xy + J1.yy*J2.yy)
end

function *(J1::HermitianJonesMatrix, J2::JonesMatrix)
    JonesMatrix(J1.xx *J2.xx + J1.xy*J2.yx, J1.xx *J2.xy + J1.xy*J2.yy,
                J1.xy'*J2.xx + J1.yy*J2.yx, J1.xy'*J2.xy + J1.yy*J2.yy)
end

function *(J1::JonesMatrix, J2::HermitianJonesMatrix)
    JonesMatrix(J1.xx*J2.xx + J1.xy*J2.xy', J1.xx*J2.xy + J1.xy*J2.yy,
                J1.yx*J2.xx + J1.yy*J2.xy', J1.yx*J2.xy + J1.yy*J2.yy)
end

\(J1::AnyJonesMatrix, J2::AnyJonesMatrix) = inv(J1)*J2
/(J1::AnyJonesMatrix, J2::AnyJonesMatrix) = J1*inv(J2)

Base.conj(J::JonesMatrix) = JonesMatrix(conj(J.xx), conj(J.xy), conj(J.yx), conj(J.yy))
Base.conj(J::DiagonalJonesMatrix) = DiagonalJonesMatrix(conj(J.xx), conj(J.yy))
Base.conj(J::HermitianJonesMatrix) = HermitianJonesMatrix(J.xx, conj(J.xy), J.yy)

Base.ctranspose(J::JonesMatrix) = JonesMatrix(conj(J.xx), conj(J.yx), conj(J.xy), conj(J.yy))
Base.ctranspose(J::DiagonalJonesMatrix) = conj(J)
Base.ctranspose(J::HermitianJonesMatrix) = J

Base.det(J::JonesMatrix) = J.xx*J.yy - J.xy*J.yx
Base.det(J::DiagonalJonesMatrix) = J.xx*J.yy
Base.det(J::HermitianJonesMatrix) = J.xx*J.yy - J.xy*J.xy'

function Base.inv(J::JonesMatrix)
    1/det(J) * JonesMatrix(J.yy, -J.xy, -J.yx, J.xx)
end

function Base.inv(J::DiagonalJonesMatrix)
    DiagonalJonesMatrix(1/J.xx, 1/J.yy)
end

function Base.inv(J::HermitianJonesMatrix)
    1/det(J) * HermitianJonesMatrix(J.yy, -J.xy, J.xx)
end

# use the Frobenius norm
Base.norm(J::JonesMatrix) = sqrt(abs2(J.xx)+abs2(J.xy)+abs2(J.yx)+abs2(J.yy))
Base.norm(J::DiagonalJonesMatrix) = sqrt(abs2(J.xx)+abs2(J.yy))
Base.norm(J::HermitianJonesMatrix) = sqrt(abs2(J.xx)+2abs2(J.xy)+abs2(J.yy))

function Base.kron(J1::JonesMatrix, J2::JonesMatrix)
    [J1.xx*J2.xx J1.xx*J2.xy J1.xy*J2.xx J1.xy*J2.xy;
     J1.xx*J2.yx J1.xx*J2.yy J1.xy*J2.yx J1.xy*J2.yy;
     J1.yx*J2.xx J1.yx*J2.xy J1.yy*J2.xx J1.yy*J2.xy;
     J1.yx*J2.yx J1.yx*J2.yy J1.yy*J2.yx J1.yy*J2.yy]
end

doc"""
    congruence_transform(J::JonesMatrix, K::HermitianJonesMatrix)

Compute the congruence transformation of $K$ with respect to $J$. Using this function instead of
explicitly computing $J*K*J'$ guarantees that the final result is exactly Hermitian.

``` math
    K \rightarrow JKJ^*
```
"""
function congruence_transform(J::JonesMatrix, K::HermitianJonesMatrix)
    HermitianJonesMatrix(abs2(J.xx)*K.xx + 2real(J.xx*J.xy'*K.xy) + abs2(J.xy)*K.yy,
                         J.xx*J.yx'*K.xx + J.xx*J.yy'*K.xy + J.xy*J.yx'*K.xy' + J.xy*J.yy'*K.yy,
                         abs2(J.yx)*K.xx + 2real(J.yx*J.yy'*K.xy) + abs2(J.yy)*K.yy)
end

function make_hermitian(J::JonesMatrix)
    HermitianJonesMatrix(real(J.xx), 0.5*(J.xy+conj(J.yx)), real(J.yy))
end

