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

#__precompile__()

module TTCal

export Nant, Nfreq, Nbase

#export JonesMatrix, DiagonalJonesMatrix, HermitianJonesMatrix
#export StokesVector, MuellerMatrix
#
#export Visibilities, Metadata, Nant, Nfreq, Nbase
#
#export Source
#export PointSource, GaussianSource, DiskSource
#export ShapeletSource, MultiSource, RFISource
#export PowerLaw, RFISpectrum
#export readsources, writesources
#export ConstantBeam, SineBeam, Memo178Beam, ZernikeBeam
#
#export genvis, subsrc!, putsrc!, getspec, fitvis
#
#export GainCalibration, gaincal
#export PolarizationCalibration, polcal
#export applycal!, corrupt!, peel!, shave!
#export PeelingSource, ShavingSource, ZestingSource, PruningSource
#
#importall Base.Operators

#using DocOpt
using ProgressMeter
using FileIO, JLD2, JSON
using NLopt # used in fitvis
using Unitful, UnitfulAstro
using CasaCore.Measures
using CasaCore.Tables
using StaticArrays

struct TTCalException <: Exception
    message :: String
end
@noinline err(message) = throw(TTCalException(message))
function Base.show(io::IO, exception::TTCalException)
    print(io, "TTCalException: ", exception.message)
end

include("special.jl")
include("rungekutta.jl")
include("jones.jl")
include("stokes.jl")
include("struct-of-arrays.jl")

#include("spectra.jl")
#include("skymodels.jl")
#include("beammodels.jl")

#include("utm.jl")
#include("ionosphere.jl")
include("metadata.jl")
include("visibilities.jl")
#include("genvis.jl")
#include("getspec.jl")
#include("subsrc.jl")
#include("fitvis.jl")
#include("calibration.jl")
#include("peel.jl")
#include("commandline.jl")

end

