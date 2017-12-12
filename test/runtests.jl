using TTCal
using Base.Test
using CasaCore.Measures
using CasaCore.Tables
using Unitful, UnitfulAstro

include("setup.jl")

srand(123)
@testset "TTCal Tests" begin
    @testset "math" begin
        include("math/special.jl")
        #include("math/rungekutta.jl")
        include("math/jones.jl")
        include("math/stokes.jl")
    end
    @testset "data" begin
    #    #include("data/metadata.jl")
    #    #include("data/visibilities.jl")
        include("data/rotate-phase-center.jl")
    end
    @testset "sky" begin
        include("sky/spectra.jl")
        #include("sky/skymodels.jl")
    end
    @testset "instrument" begin
        include("instrument/beams.jl")
        include("instrument/genvis.jl")
    end
    #include("utm.jl")
    #include("ionosphere.jl")
    include("getspec.jl")
    #include("subsrc.jl")
    include("fitvis.jl")
    include("calibration.jl")
    #include("peel.jl")
end

