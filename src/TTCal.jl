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

__precompile__()

module TTCal

export PointSource
export readsources, writesources

export genvis
export getspec
export fitvis
export subsrc!

export AmplitudeCalibration, ampcal
export GainCalibration, gaincal
export PolarizationCalibration, polcal
export applycal!, corrupt!
export peel!

importall Base.Operators
import Base: zero, one, rand, det, inv, norm

using JSON
using JLD
using CasaCore.Measures
using CasaCore.Tables
using MeasurementSets

const c = 2.99792e+8
include("rungekutta.jl")
include("jones.jl")
include("sourcemodel.jl")
include("fringepattern.jl")
include("genvis.jl")
include("getspec.jl")
include("fitvis.jl")
include("subsrc.jl")
include("calibration.jl")
include("peel.jl")
include("beammodel.jl")
include("utm.jl")

function run_gaincal(args)
    ms = Table(ascii(args["--input"]))
    sources = readsources(args["--sources"])
    maxiter = haskey(args,"--maxiter")? args["--maxiter"] : 20
    tolerance = haskey(args,"--tolerance")? args["--tolerance"] : 1e-4
    force_imaging_columns = haskey(args,"--force-imaging")
    cal = gaincal(ms,sources,
                  maxiter=maxiter,
                  tolerance=tolerance,
                  force_imaging_columns=force_imaging_columns)
    write(args["--output"],cal)
    cal
end

function run_polcal(args)
    ms = Table(ascii(args["--input"]))
    sources = readsources(args["--sources"])
    maxiter = haskey(args,"--maxiter")? args["--maxiter"] : 20
    tolerance = haskey(args,"--tolerance")? args["--tolerance"] : 1e-4
    force_imaging_columns = haskey(args,"--force-imaging")
    cal = polcal(ms,sources,
                 maxiter=maxiter,
                 tolerance=tolerance,
                 force_imaging_columns=force_imaging_columns)
    write(args["--output"],cal)
    cal
end

function run_peel(args)
    ms = Table(ascii(args["--input"]))
    sources = readsources(args["--sources"])
    minuvw = haskey(args,"--minuvw")? args["--minuvw"] : 15.0
    peel!(GainCalibration,ms,sources,minuvw=minuvw)
end

function run_applycal(args)
    cal = read(args["--calibration"])
    force_imaging_columns = haskey(args,"--force-imaging")
    apply_to_corrected = haskey(args,"--corrected")
    for input in args["--input"]
        ms = Table(ascii(input))
        applycal!(ms,cal,
                  force_imaging_columns=force_imaging_columns,
                  apply_to_corrected=apply_to_corrected)
    end
    cal
end

end

