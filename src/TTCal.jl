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

export Source
export genvis
export fitvis
export subsrc!

export clearflags!,flagdata!
export bandpass
export polcal
export applycal!

using JSON, NPZ
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

include("flagdata.jl")
include("bandpass.jl")
include("polcal.jl")
include("applycal.jl")

function run_flagdata(args)
    for input in args["--input"]
        ms = Table(ascii(input))
        flagdata!(ms,args["--antennas"])
    end
    nothing
end

function run_bandpass(args)
    ms = Table(ascii(args["--input"]))
    sources = readsources(args["--sources"])
    maxiter = haskey(args,"--maxiter")? args["--maxiter"] : 20
    tol = haskey(args,"--tolerance")? args["--tolerance"] : 1e-4
    criteria = StoppingCriteria(maxiter,tol)
    gains = bandpass(ms,sources,criteria)
    npzwrite(args["--output"],gains)
    nothing
end

function run_applycal(args)
    gains = npzread(args["--calibration"])
    for input in args["--input"]
        ms = Table(ascii(input))
        applycal!(ms,gains)
    end
    nothing
end

end

