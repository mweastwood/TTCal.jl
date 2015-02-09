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
        ms = Table(input)
        flagdata!(ms,args["--antennas"])
    end
    nothing
end

function run_bandpass(args)
    ms = Table(args["--input"])
    sources = readsources(args["--sources"])
    criteria = StoppingCriteria(20,1e-4)
    gains = bandpass(ms,sources,criteria)
    npzwrite(args["--output"],gains)
    nothing
end

function run_applycal(args)
    gains = npzread(args["--calibration"])
    for input in args["--input"]
        ms = Table(input)
        applycal!(ms,gains)
    end
    nothing
end

end

