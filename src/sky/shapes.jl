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

abstract type AbstractShape{Spectrum} end

struct Point{Spectrum} <: AbstractShape{Spectrum}
    direction :: Direction
    spectrum  :: Spectrum
end

struct Gaussian{Spectrum} <: AbstractShape{Spectrum}
    direction :: Direction
    spectrum  :: Spectrum
    major_fwhm     :: Float64 # FWHM along the major axis (radians)
    minor_fwhm     :: Float64 # FWHM along the minor axis (radians)
    position_angle :: Float64 # (radians)
end

struct Disk{Spectrum} <: AbstractShape{Spectrum}
    direction :: Direction
    spectrum  :: Spectrum
    radius    :: Float64 # radius of the disk (radians)
end

struct Shapelet{Spectrum} <: AbstractShape{Spectrum}
    direction :: Direction
    spectrum  :: Spectrum
    scale :: Float64         # scale factor of the shapelets (radians)
    coeff :: Vector{Float64} # list of shapelet coefficients
end

#doc"""
#    ShapeletSource <: Source
#
#An astronomical source composed of shapelets.
#
#Shapelets are eigenfunctions of the fourier transform operator.  They form an orthonormal basis for
#real functions mapping $\mathbb R^2 \mapsto \mathbb R$.
#"""
#"""
#    RFISource <: Source
#
#A terrestrial source of RFI. These sources are assumed to be spectrally unsmooth and in the near
#field of the interferometer.
#"""
#type RFISource <: Source
#    name :: String
#    position :: Position
#    spectrum :: RFISpectrum
#end




function isabovehorizon(frame::ReferenceFrame, direction::Direction; threshold=0)
    azel = measure(frame, direction, dir"AZEL")
    el = latitude(azel)
    el > threshold
end

function isabovehorizon(frame::ReferenceFrame, shape::AbstractShape; threshold=0)
    isabovehorizon(frame, shape.direction, threshold)
end

function isrising(frame::ReferenceFrame, shape::AbstractShape)
    azel = measure(frame, shape.direction, dir"AZEL")
    az = longitude(azel) # ∈ (-π/2, π/2)
    az > 0
end

function Base.:*(number::Real, shape::Point)
    Point(shape.direction, number*shape.spectrum)
end

function Base.:*(number::Real, shape::Gaussian)
    Gaussian(shape.direction, number*shape.spectrum,
             shape.major_fwhm, shape.minor_fwhm, shape.position_angle)
end

Base.:*(shape::AbstractShape, number::Real) = number*shape

