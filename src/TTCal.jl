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

end

