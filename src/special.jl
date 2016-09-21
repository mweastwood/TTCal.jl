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
    hermite(n, x)

Compute the value of the $n$th Hermite polynomial at $x$.
"""
function hermite(n, x)
    if n == 0
        return one(typeof(x))
    elseif n == 1
        return 2x
    else
        # Note that the naive thing to do here is:
        #
        #     return 2x*hermite(n-1, x) - 2*(n-1)*hermite(n-2, x)
        #
        # This is sub-optimal because you end up calculating `hermite(n-2, x)`
        # two times, and `hermite(n-3, x)` three times, etc. So instead of
        # using recursion we will just write out the for-loop because you end
        # up doing less work (especially for large-n).
        local output
        minus_two = hermite(0, x)
        minus_one = hermite(1, x)
        for n′ = 2:n
            output = 2x*minus_one - 2*(n′-1)*minus_two
            minus_two = minus_one
            minus_one = output
        end
        return output
    end
end

doc"""
    zernike(n, m, ρ, θ)

Compute the value of the Zernike polynomial $Z_{n,m}$ at the polar coordinates $(ρ,θ)$.
"""
function zernike(n, m, ρ, θ)
    zernike_radial_part(n, abs(m), ρ) * zernike_azimuthal_part(m, θ)
end

function zernike_radial_part(n, m, ρ)
    R0 = ρ^m
    n == m && return R0
    R2 = ((m+2)*ρ^2 - (m+1))*R0
    for n′ = m+4:2:n
        recurrence_relation = ((2*(n′-1)*(2n′*(n′-2)*ρ^2-m^2-n′*(n′-2))*R2 - n′*(n′+m-2)*(n′-m-2)*R0)
                                    / ((n′+m)*(n′-m)*(n′-2)))
        R0 = R2
        R2 = recurrence_relation
    end
    R2
end

function zernike_azimuthal_part(m, θ)
    if m == 0
        return 1.0
    elseif m > 0
        return cos(m*θ)
    else
        return sin(-m*θ)
    end
end

