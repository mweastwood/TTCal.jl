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

"""
This type represents a 2x2 complex Jones matrix.

    ⌈xx xy⌉
    ⌊yx yy⌋
"""
immutable JonesMatrix
    xx::Complex64
    xy::Complex64
    yx::Complex64
    yy::Complex64
end

JonesMatrix() = one(JonesMatrix)
zero(::Type{JonesMatrix}) = JonesMatrix(0,0,0,0)
one(::Type{JonesMatrix}) = JonesMatrix(1,0,0,1) # the identity matrix

for op in (:+,:-)
    @eval function $op(J1::JonesMatrix,J2::JonesMatrix)
        JonesMatrix($op(J1.xx,J2.xx),
                    $op(J1.xy,J2.xy),
                    $op(J1.yx,J2.yx),
                    $op(J1.yy,J2.yy))
    end
end

function *(a::Number,J::JonesMatrix)
    JonesMatrix(a*J.xx,a*J.xy,a*J.yx,a*J.yy)
end
*(J::JonesMatrix,a::Number) = *(a,J)

@inline function *(J1::JonesMatrix,J2::JonesMatrix)
    JonesMatrix(J1.xx*J2.xx + J1.xy*J2.yx,
                J1.xx*J2.xy + J1.xy*J2.yy,
                J1.yx*J2.xx + J1.yy*J2.yx,
                J1.yx*J2.xy + J1.yy*J2.yy)
end

function \(J1::JonesMatrix,J2::JonesMatrix)
    inv(J1)*J2
end

function ctranspose(J::JonesMatrix)
    JonesMatrix(conj(J.xx),conj(J.yx),conj(J.xy),conj(J.yy))
end

function det(J::JonesMatrix)
    J.xx*J.yy - J.xy*J.yx
end

function inv(J::JonesMatrix)
    d = det(J)
    JonesMatrix(J.yy/d,-J.xy/d,-J.yx/d,J.xx/d)
end

function norm(J::JonesMatrix)
    # The Frobenius norm
    sqrt(abs2(J.xx)+abs2(J.xy)+abs2(J.yx)+abs2(J.yy))
end

