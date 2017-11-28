using TTCal
using Base.Test
using CasaCore.Measures
using CasaCore.Tables

include("setup.jl")

srand(123)
@testset "TTCal Tests" begin
    @testset "math" begin
        include("math/special.jl")
        #include("math/rungekutta.jl")
        include("math/jones.jl")
        #include("math/stokes.jl")
    end
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

