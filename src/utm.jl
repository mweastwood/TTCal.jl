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

immutable Hemisphere{N} end
const north = Hemisphere{true}()
const south = Hemisphere{false}()

"""
    WGS84

A module containing constants that define the WGS 84 spatial
reference system.
"""
module WGS84
    const a = 6378.137 # equatorial radius (km)
    const f = 1/298.257223563 # flattening
end

"""
    latlong_to_utm(zone, latitude, longitude)

Convert the latitude and longitude (both in degrees) to
an easting and a northing (both in meters).

This function is entirely based on the [Universal Transverse Mercator coordinate system]
(https://en.wikipedia.org/w/index.php?title=Universal_Transverse_Mercator_coordinate_system&oldid=683193579#Simplified_formulas)
Wikipedia page.
"""
function latlong_to_utm(zone, latitude, longitude)
    k0 = 0.9996
    N0 = latitude ≥ 0? 0.0 : 10000.0
    E0 = 500

    n = WGS84.f/(2-WGS84.f)
    A = WGS84.a/(1+n) * (1 + n^2/4 + n^4/64)

    α = [1/2  -2/3   5/16;
           0 13/48   -3/5;
           0     0 61/240] * [n;n^2;n^3]

    λ0 = (zone*6 - 183)*π/180
    λ = longitude*π/180
    ϕ =  latitude*π/180

    t = sinh(atanh(sin(ϕ)) - 2*sqrt(n)/(1+n)*atanh(2*sqrt(n)/(1+n)*sin(ϕ)))
    ξ′ = atan(t/cos(λ-λ0))
    η′ = atanh(sin(λ-λ0)/sqrt(1+t^2))

    σ = 1.0
    τ = 0.0
    for j = 1:3
        σ += 2*j*α[j]*cos(2*j*ξ′)*cosh(2*j*η′)
        τ += 2*j*α[j]*sin(2*j*ξ′)*sinh(2*j*η′)
    end

    easting  = E0 + k0*A*η′
    northing = N0 + k0*A*ξ′
    for j = 1:3
        easting  += k0*A*α[j]*cos(2*j*ξ′)*sinh(2*j*η′)
        northing += k0*A*α[j]*sin(2*j*ξ′)*cosh(2*j*η′)
    end
    easting  *= 1e3
    northing *= 1e3

    easting,northing
end

"""
    utm_to_latlong(zone, easting, northing, hemisphere = north)

Convert the easting and northing (both in meters) to
a latitude and longitude (both in degrees).

This function is entirely based on the [Universal Transverse Mercator coordinate system]
(https://en.wikipedia.org/w/index.php?title=Universal_Transverse_Mercator_coordinate_system&oldid=683193579#Simplified_formulas)
Wikipedia page.
"""
function utm_to_latlong{N}(zone, easting, northing,
                           hemisphere::Hemisphere{N} = north)
    k0 = 0.9996
    N0 = N? 0.0 : 10000.0
    E0 = 500

    n = WGS84.f/(2-WGS84.f)
    A = WGS84.a/(1+n) * (1 + n^2/4 + n^4/64)

    β = [1/2  -2/3  37/96;
           0  1/48   1/15;
           0     0 17/480] * [n;n^2;n^3]
    δ = [2 -2/3    -2;
         0  7/3  -8/5;
         0    0 56/15] * [n;n^2;n^3]
    
    ξ = (northing/1e3-N0)/(k0*A)
    η = ( easting/1e3-E0)/(k0*A)
    ξ′ = ξ
    η′ = η
    σ′ = 1.0
    τ′ = 0.0
    for j = 1:3
        ξ′ -= β[j]*sin(2*j*ξ)*cosh(2*j*η)
        η′ -= β[j]*cos(2*j*ξ)*sinh(2*j*η)
        σ′ -= 2*j*β[j]*cos(2*j*ξ)*cosh(2*j*η)
        τ′ += 2*j*β[j]*sin(2*j*ξ)*sinh(2*j*η)
    end

    λ0 = zone*6 - 183
    longitude = λ0 + atan(sinh(η′)/cos(ξ′))*180/π

    χ = asin(sin(ξ′)/cosh(η′))
    latitude = χ
    for j = 1:3
        latitude += δ[j]*sin(2*j*χ)
    end
    latitude *= 180/π

    latitude,longitude
end

"""
    azel_to_itrf(direction, location)

The conversion from azimuth and elevation to an ITRF direction is
straightforward. CasaCore's implementation is slow and requires
unnecessary information (the current time).
"""
function azel_to_itrf(direction, location)
    direction.sys == dir"AZEL" || error("Direction must be in the AZEL coordinate system.")
    location.sys  == pos"ITRF" || error("Location must be in the ITRF coordinate system.")
    z     = [0, 0, 1]
    loc   = [location.x, location.y, location.z]
    up    = loc / norm(loc)
    north = z - dot(z, up)*up
    north = north / norm(north)
    east  = cross(north, up)
    dir   = direction.x * north + direction.y * east + direction.z * up
    Direction(dir"ITRF", dir[1], dir[2], dir[3])
end

