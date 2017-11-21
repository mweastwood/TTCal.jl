@testset "skymodels.jl" begin
    frame   = ReferenceFrame()
    sources = [PointSource("S1",Direction(dir"AZEL","0.0deg","+45.0deg"),
                                PowerLaw(1.0,2.0,3.0,4.0,10e6,[0.0])),
               PointSource("S2",Direction(dir"AZEL","0.0deg","-45.0deg"),
                                PowerLaw(1.0,2.0,3.0,4.0,10e6,[0.0]))]
    @test TTCal.isabovehorizon(frame,sources[1]) == true
    @test TTCal.isabovehorizon(frame,sources[2]) == false
    sources′ = TTCal.abovehorizon(frame,sources)
    @test length(sources′) == 1
    @test sources′[1].name == "S1"

    sources = readsources("sources.json")
    cas_a = sources[1]
    cyg_a = sources[2]
    @test cas_a == PointSource("Cas A",Direction(dir"J2000","23h23m24s","58d48m54s"),
                                       PowerLaw(555904.26,0,0,0,1e6,[-0.770]))
    @test cyg_a == PointSource("Cyg A",Direction(dir"J2000","19h59m28.35663s","+40d44m02.0970s"),
                                       PowerLaw(49545.02,0,0,0,1e6,[+0.085,-0.178]))

    sources1 = [PointSource("S1",Direction(dir"J2000","1h","0d"),
                                 PowerLaw(1.0,2.0,3.0,4.0,10e6,[0.0])),
                PointSource("S2",Direction(dir"J2000","2h","0d"),
                                 PowerLaw(1.0,2.0,3.0,4.0,10e6,[0.0])),
                PointSource("S3",Direction(dir"J2000","3h","0d"),
                                 PowerLaw(1.0,2.0,3.0,4.0,10e6,[0.0])),
                GaussianSource("gaussian", Direction(dir"J2000", "4h", "10d"),
                                           PowerLaw(4, 3, 2, 1, 10e6, [1.0]), 0, 0, 0)]
    filename = tempname()*".json"
    writesources(filename,sources1)
    sources2 = readsources(filename)
    @test sources1 == sources2

    multisources1 = [MultiSource("M", [PointSource("S1",Direction(dir"J2000","1h","0d"),
                                                   PowerLaw(1.0,2.0,3.0,4.0,10e6,[0.0])),
                                       PointSource("S2",Direction(dir"J2000","2h","0d"),
                                                   PowerLaw(1.0,2.0,3.0,4.0,10e6,[0.0])),
                                       PointSource("S3",Direction(dir"J2000","3h","0d"),
                                                   PowerLaw(1.0,2.0,3.0,4.0,10e6,[0.0]))])]

    filename = tempname()*".json"
    writesources(filename,multisources1)
    multisources2 = readsources(filename)
    @test multisources1 == multisources2

    @testset "rfi sources" begin
        # Dear LADWP,
        #
        # Please fix.
        #
        # Love,
        # Michael
        pos = Position(pos"WGS84", 1226.709meters, -118.31478degrees, 37.14540degrees)
        ν = collect(linspace(30e6, 60e6, 5))
        stokes = [rand(StokesVector) for idx = 1:length(ν)]
        pls_fix_this_ladwp = RFISource("Power Lines", pos, RFISpectrum(ν, stokes))
        filename = tempname()*".json"
        writesources(filename, [pls_fix_this_ladwp])
        ladwp_pls = readsources(filename)[1]
        @test pls_fix_this_ladwp.name == ladwp_pls.name
        @test pls_fix_this_ladwp.position ≈ ladwp_pls.position
        @test pls_fix_this_ladwp.spectrum.channels == ladwp_pls.spectrum.channels
        @test pls_fix_this_ladwp.spectrum.stokes == ladwp_pls.spectrum.stokes

        ν = collect(linspace(30e6, 60e6, 5))
        stokes = [rand(StokesVector) for idx = 1:length(ν)+1]
        @test_throws ErrorException RFISource("too many stokes vectors", pos, RFISpectrum(ν, stokes))
    end
end

