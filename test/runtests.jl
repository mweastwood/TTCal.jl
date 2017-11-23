using TTCal
using Base.Test
using CasaCore.Measures
using CasaCore.Tables
#using NPZ

include("setup.jl")

srand(123)
@testset "TTCal Tests" begin
    include("special.jl")
    #include("rungekutta.jl")
    include("jones.jl")
    #include("stokes.jl")
    #include("spectra.jl")
    #include("skymodels.jl")
    #include("beammodels.jl")
    #include("utm.jl")
    #include("ionosphere.jl")
    #include("metadata.jl")
    #include("visibilities.jl")
    #include("genvis.jl")
    #include("getspec.jl")
    #include("subsrc.jl")
    #include("fitvis.jl")
    #include("calibration.jl")
    #include("peel.jl")
end

