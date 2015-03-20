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

include("bandpass.jl")
include("polcal.jl")
include("applycal.jl")

function run_bandpass(args)
    ms = Table(ascii(args["--input"]))
    sources = readsources(args["--sources"])
    maxiter = haskey(args,"--maxiter")? args["--maxiter"] : 20
    tol = haskey(args,"--tolerance")? args["--tolerance"] : 1e-4
    criteria = StoppingCriteria(maxiter,tol)
    gains,gain_flags = bandpass(ms,sources,criteria)

    # Write the gains to a file
    Nant,Npol,Nchan = size(gains)
    open(args["--output"],"w") do f
        write(f,'B')
        write(f,Int32(Nant))
        write(f,Int32(Nchan))
        write(f,permutedims(gain_flags,(2,3,1)))
        write(f,permutedims(gains,(2,3,1)))
    end

    gains, gain_flags
end

function run_polcal(args)
    ms = Table(ascii(args["--input"]))
    sources = readsources(args["--sources"])
    maxiter = haskey(args,"--maxiter")? args["--maxiter"] : 20
    tol = haskey(args,"--tolerance")? args["--tolerance"] : 1e-4
    criteria = StoppingCriteria(maxiter,tol)
    gains,gain_flags = polcal(ms,sources,criteria)

    # Write the gains to a file
    Npol1,Npol2,Nant,Nchan = size(gains)
    open(args["--output"],"w") do f
        write(f,'J')
        write(f,Int32(Nant))
        write(f,Int32(Nchan))
        write(f,permutedims(gain_flags,(2,1)))
        write(f,permutedims(gains,(1,2,4,3)))
    end

    gains, gain_flags
end

function run_applycal(args)
    local T, Nant, Nchan, gain_flags, gains

    # Read in the gains
    open(args["--calibration"],"r") do f
        T = read(f,Char)
        Nant = Int(read(f,Int32))
        Nchan = Int(read(f,Int32))
        if T == 'B'
            gain_flags = permutedims(read(f,Bool,(2,Nchan,Nant)),(3,1,2))
            gains = permutedims(read(f,Complex64,(2,Nchan,Nant)),(3,1,2))
        elseif T == 'J'
            gain_flags = permutedims(read(f,Bool,(Nchan,Nant)),(2,1))
            gains = permutedims(read(f,Complex64,(2,2,Nchan,Nant)),(1,2,4,3))
        else
            error("Unknown calibration type.")
        end
    end

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

