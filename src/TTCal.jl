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

module TTCal

export PointSource
export readsources, writesources

export genvis
export getspec
export fitvis
export subsrc!

export bandpass
export polcal
export applycal!

using JSON
using SIUnits
using CasaCore.Measures
using CasaCore.Tables

include("units.jl")
include("ms.jl")
include("rungekutta.jl")

include("sourcemodel.jl")
include("fringepattern.jl")
include("genvis.jl")
include("getspec.jl")
include("fitvis.jl")
include("subsrc.jl")

include("io.jl")
include("bandpass.jl")
include("polcal.jl")
include("applycal.jl")

function run_bandpass(args)
    ms = Table(ascii(args["--input"]))
    sources = haskey(args,"--sources")? readsources(args["--sources"]) : Source[]
    maxiter = haskey(args,"--maxiter")? args["--maxiter"] : 20
    tol = haskey(args,"--tolerance")? args["--tolerance"] : 1e-4
    criteria = StoppingCriteria(maxiter,tol)
    force_imaging_columns = haskey(args,"--force-imaging")
    model_already_present = !haskey(args,"--sources")
    gains,gain_flags = bandpass(ms,sources,criteria,
                                force_imaging_columns=force_imaging_columns,
                                model_already_present=model_already_present)
    write_gains(args["--output"],gains,gain_flags)
    gains, gain_flags
end

function run_polcal(args)
    ms = Table(ascii(args["--input"]))
    sources = haskey(args,"--sources")? readsources(args["--sources"]) : Source[]
    maxiter = haskey(args,"--maxiter")? args["--maxiter"] : 20
    tol = haskey(args,"--tolerance")? args["--tolerance"] : 1e-4
    criteria = StoppingCriteria(maxiter,tol)
    force_imaging_columns = haskey(args,"--force-imaging")
    model_already_present = !haskey(args,"--sources")
    gains,gain_flags = polcal(ms,sources,criteria,
                              force_imaging_columns=force_imaging_columns,
                              model_already_present=model_already_present)
    write_gains(args["--output"],gains,gain_flags)
    gains, gain_flags
end

function run_applycal(args)
    gains,gain_flags = read_gains(args["--calibration"])
    force_imaging_columns = haskey(args,"--force-imaging")
    apply_to_corrected = haskey(args,"--corrected")
    for input in args["--input"]
        ms = Table(ascii(input))
        applycal!(ms,gains,gain_flags,
                  force_imaging_columns=force_imaging_columns,
                  apply_to_corrected=apply_to_corrected)
    end
    nothing
end

end

