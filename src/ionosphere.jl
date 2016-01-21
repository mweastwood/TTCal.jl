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
    appleton_hartree(frequency, plasma_frequency)

Compute the index of refraction of the ionosphere.

We currently assume a cold, collisionless, and unmagnetized plasma.
"""
function appleton_hartree(frequency, plasma_frequency)
    ω  = 2π*frequency
    ωp = 2π*plasma_frequency
    X  = ωp^2 / ω^2
    sqrt(1 - X)
end

const Rearth = 6300e3 # radius of the Earth
# The F layer of the ionosphere is roughly 200 km to 400 km.
const Flayer = (200e3, 400e3)

"""
    ray_trace(elevation, frequency, plasma_frequency)

Trace a ray through the ionosphere to determine the apparent
position of a source.
"""
function ray_trace(elevation, frequency, plasma_frequency)
    index_of_refraction = appleton_hartree(frequency, plasma_frequency)
    # Start a ray at the surface of the earth
    θ = elevation
    m = tan(θ)
    b = Rearth
    # Refract off the inner edge of the ionosphere
    x,y = find_intersection_with_circle(m, b, Rearth+Flayer[1])
    angle_of_incidence  = atan2(y, x) - θ
    angle_of_refraction = asin(sin(angle_of_incidence)/index_of_refraction)
    θ = atan2(y, x) - angle_of_refraction
    m = tan(θ)
    b = y - m*x
    # Refract off the outer edge of the ionosphere
    x,y = find_intersection_with_circle(m, b, Rearth+Flayer[2])
    angle_of_incidence  = atan2(y, x) - θ
    angle_of_refraction = asin(sin(angle_of_incidence)*index_of_refraction)
    θ = atan2(y, x) - angle_of_refraction
end

doc"""
    find_intersection_with_circle(m, b, R)

Find the coordinates $(x,y)$ where the line $y = mx+b$
intersects the circle $x^2+y^2 = R^2$.
"""
function find_intersection_with_circle(m, b, R)
    A = 1 + m^2
    B = 2m*b
    C = b^2 - R^2
    x = (-B + sqrt(B^2 - 4A*C))/(2A)
    y = m*x + b
    x,y
end

