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
    gains,gain_flags = bandpass(ms,sources,criteria)

    # Write the gains to a CASA .bcal table
    Nant,Npol,Nchan = size(gains)
    bcal = Table(ascii(args["--output"]))
    Nrows = numrows(bcal)
    if Nrows < Nant
        Tables.addRows!(bcal,Nant-Nrows)
    elseif Nrows > Nant
        Tables.removeRows!(bcal,[Nant+1:Nrows...])
    end
    bcal["ANTENNA1"] = Cint[1:Nant...]
    bcal["CPARAM"] = permutedims(gains,(2,3,1))
    bcal["FLAG"] = permutedims(gain_flags,(2,3,1))
    gains
end

function run_applycal(args)
    bcal = Table(ascii(args["--calibration"]))
    gains = permutedims(bcal["CPARAM"],(3,1,2))
    gain_flags = permutedims(bcal["FLAG"],(3,1,2))
    for input in args["--input"]
        ms = Table(ascii(input))
        applycal!(ms,gains,gain_flags)
    end
    nothing
end

end

