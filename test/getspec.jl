function test_getspec_direction(metadata, source, direction)
    data = genvis(metadata, TTCal.ConstantBeam(), source)
    stokes_flux   = TTCal.total_flux.(source, metadata.frequencies)
    measured_flux = getspec(data, direction)
    @test stokes_flux ≈ measured_flux
end

function test_getspec_source(metadata, source)
    data = genvis(metadata, TTCal.ConstantBeam(), source)
    stokes_flux   = TTCal.total_flux.(source, metadata.frequencies)
    measured_flux = getspec(data, source)
    @test stokes_flux ≈ measured_flux
end

@testset "getspec.jl" begin
    metadata = TTCal.Metadata([50.0u"MHz", 74.3u"MHz", 93.1u"MHz"],
                              [Epoch(epoch"UTC", (2017-1858)*u"yr")],
                              [Position(pos"ITRF", 100randn(), 100randn(), 100randn()) for n = 1:3],
                              [Direction(dir"AZEL", 0u"°", 90u"°")])


    @testset "point sources" begin
        # unpolarized
        source = TTCal.Source("FRB",
                              TTCal.Point(Direction(dir"AZEL", 10u"°", 30u"°"),
                                          TTCal.PowerLaw(rand(), 0, 0, 0, 100u"MHz", [-rand()])))
        test_getspec_direction(metadata, source, Direction(dir"AZEL", 10u"°", 30u"°"))
        test_getspec_source(metadata, source)

        # polarized
        stokes = rand(TTCal.StokesVector)
        source = TTCal.Source("FRB",
                              TTCal.Point(Direction(dir"AZEL", 10u"°", 30u"°"),
                                          TTCal.PowerLaw(stokes, 100u"MHz", [-rand()])))
        test_getspec_direction(metadata, source, Direction(dir"AZEL", 10u"°", 30u"°"))
        test_getspec_source(metadata, source)
    end

    @testset "gaussian sources" begin
        stokes = rand(TTCal.StokesVector)
        source = TTCal.Source("test",
                              TTCal.Gaussian(Direction(dir"AZEL", 10u"°", 30u"°"),
                                             TTCal.PowerLaw(stokes, 10u"MHz", [-rand()]),
                                             deg2rad(500/3600), deg2rad(300/3600), deg2rad(50)))
        test_getspec_source(metadata, source)
    end

    @testset "multi-component sources" begin
        # this test is pretty gnarly
        stokes = rand(TTCal.StokesVector)
        point = TTCal.Point(Direction(dir"AZEL", 10u"°", 30u"°"),
                            TTCal.PowerLaw(stokes, 10u"MHz", [-rand()]))
        stokes = rand(TTCal.StokesVector)
        gaussian = TTCal.Gaussian(Direction(dir"AZEL", 10u"°", 30u"°"),
                                  TTCal.PowerLaw(stokes, 10u"MHz", [-rand()]),
                                  deg2rad(500/3600), deg2rad(300/3600), deg2rad(50))
        source = TTCal.Source("multi", [point, gaussian])
        test_getspec_source(metadata, source)
    end
end

