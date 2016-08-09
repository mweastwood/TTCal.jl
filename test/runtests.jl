using TTCal
if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end
using CasaCore.Measures
using CasaCore.Tables
using NPZ

include("setup.jl")

srand(123)
@testset "TTCal Tests" begin
    include("special.jl")
    include("stokes.jl")
    include("sourcemodel.jl")
    include("beammodel.jl")
    include("utm.jl")
    include("ionosphere.jl")
    include("measurementsets.jl")
    include("genvis.jl")
    include("getspec.jl")
    include("subsrc.jl")
    include("fitvis.jl")
    include("calibration.jl")
    include("peel.jl")
end

