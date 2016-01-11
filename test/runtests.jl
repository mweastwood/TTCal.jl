using TTCal
if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end
using CasaCore.Measures
using CasaCore.Tables

include("setup.jl")

srand(123)
#include("measurementsets.jl")
#include("stokes.jl")
#include("sourcemodel.jl")
include("beammodel.jl")
#include("fringepattern.jl")
#include("getspec.jl")
#include("genvis.jl")
#include("subsrc.jl")
#include("fitvis.jl")
#include("calibration.jl")
#include("peel.jl")
#include("utm.jl")

