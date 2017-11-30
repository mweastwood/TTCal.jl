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

export Nfreq, Ntime, Nant, Nbase
export readsky, genvis, calibrate, applycal!

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
    print(io, "TTCal Error: ", exception.message)
end

abstract type Polarization end
struct XX   <: Polarization end
struct XY   <: Polarization end
struct YX   <: Polarization end
struct YY   <: Polarization end
struct Dual <: Polarization end
struct Full <: Polarization end
const Single = Union{XX, XY, YX, YY}

include("math/special.jl")
include("math/rungekutta.jl")
include("math/jones.jl")
include("math/stokes.jl")
include("math/struct-of-arrays.jl")

include("data/metadata.jl")
include("data/visibilities.jl")

include("sky/spectra.jl")
include("sky/shapes.jl")
include("sky/sources.jl")
include("sky/models.jl")

include("instrument/beams.jl")
include("instrument/genvis.jl")

#include("sky/models.jl")
#include("sky/skymodels.jl")


#include("beammodels.jl")

#include("utm.jl")
#include("ionosphere.jl")

#include("getspec.jl")
#include("subsrc.jl")
#include("fitvis.jl")
include("calibration.jl")
#include("peel.jl")
#include("commandline.jl")

end

