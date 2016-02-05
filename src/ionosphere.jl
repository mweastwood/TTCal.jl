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

type Ionosphere
    TEC :: Float64
end

function plasma_frequency(number_density)
    _2π = 1/(2π)
    charge_of_electron = 1.6021766e-19 # C
    mass_of_electron   = 9.10938e-31 # kg
    electric_constant  = 8.854e-12 # F/m (SI units)
    _2π * charge_of_electron * sqrt(number_density / (mass_of_electron * electric_constant))
end

function tec_to_density(tec)
    # note that one TEC unit is 10^16 / m^2
    number_density = 1e16 * tec / (Flayer[2] - Flayer[1])
end

"""
    appleton_hartree(frequency, plasma_frequency)

Compute the index of refraction of the ionosphere assuming
a cold, collisionless, and unmagnetized plasma.
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

Trace a ray outwards through the ionosphere to compute
the outgoing elevation angle.
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
    θ::Float64
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

"""
    refract(elevation, frequency, plasma_frequency)

Compute the apparent elevation of a source above the ionosphere
given its original elevation.
"""
function refract(elevation, frequency, plasma_frequency)
    # use the secant method to invert the ray_trace function
    tol = deg2rad(0.1/3600) # 0.1 arcseconds
    iter = 0; maxiter = 10
    x1 = elevation
    y1 = ray_trace(x1, frequency, plasma_frequency) - elevation
    x2 = elevation + deg2rad(0.1)
    y2 = ray_trace(x2, frequency, plasma_frequency) - elevation
    while abs(x2 - x1) > tol && iter < maxiter
        xnew = min((x1*y2 - x2*y1) / (y2 - y1), π/2)
        x1 = x2
        y1 = y2
        x2 = xnew
        y2 = ray_trace(x2, frequency, plasma_frequency) - elevation
        iter += 1
    end
    x2::Float64
end

function refract(direction::Direction, frequency, plasma_frequency)
    direction.sys == dir"AZEL" || error("Direction must be in the AZEL coordinate system.")
    azimuth   = longitude(direction)
    elevation =  latitude(direction)
    elevation = refract(elevation, frequency, plasma_frequency)
    Direction(dir"AZEL", azimuth*radians, elevation*radians)
end

