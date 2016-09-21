function test_getspec_direction(meta, source)
    data = genvis(meta, source)
    dir = source.direction
    stokes_flux = source.spectrum.(meta.channels)
    measured_flux = getspec(data, meta, dir)
    for β = 1:Nfreq(meta)
        linear_flux = HermitianJonesMatrix(stokes_flux[β])
        @test linear_flux.xx ≈ measured_flux[β].xx
        @test linear_flux.xy ≈ measured_flux[β].xy
        @test linear_flux.yy ≈ measured_flux[β].yy
    end
end

function test_getspec_source(meta, source)
    data = genvis(meta, source)
    stokes_flux = source.spectrum.(meta.channels)
    measured_flux = getspec(data, meta, source)
    for β = 1:Nfreq(meta)
        linear_flux = HermitianJonesMatrix(stokes_flux[β])
        @test linear_flux.xx ≈ measured_flux[β].xx
        @test linear_flux.xy ≈ measured_flux[β].xy
        @test linear_flux.yy ≈ measured_flux[β].yy
    end
end

function test_getspec_multi(meta, source)
    data = genvis(meta, source)
    linear_flux = [TTCal.get_total_flux(source, ν) for ν in meta.channels]
    measured_flux = getspec(data, meta, source)
    for β = 1:Nfreq(meta)
        @test linear_flux[β].xx ≈ measured_flux[β].xx
        @test linear_flux[β].xy ≈ measured_flux[β].xy
        @test linear_flux[β].yy ≈ measured_flux[β].yy
    end
end

@testset "getspec.jl" begin
    Nant = 10
    Nfreq = 2
    name, ms = createms(Nant, Nfreq)
    meta = Metadata(ms)
    beam = ConstantBeam()

    @testset "point sources" begin
        # unpolarized
        source = PointSource("FRB", Direction(dir"AZEL", 10degrees, 30degrees),
                                    PowerLaw(rand(), 0, 0, 0, 10e6, [-rand()]))
        test_getspec_direction(meta, source)
        test_getspec_source(meta, source)

        # polarized
        source = PointSource("FRB", Direction(dir"AZEL", 10degrees, 30degrees),
                                    PowerLaw(rand(StokesVector), 10e6, [-rand()]))
        test_getspec_direction(meta, source)
        test_getspec_source(meta, source)
    end

    @testset "gaussian sources" begin
        source = GaussianSource("test", Direction(dir"AZEL", 10degrees, 30degrees),
                                        PowerLaw(rand(StokesVector), 10e6, [-rand()]),
                                        deg2rad(500/3600), deg2rad(300/3600), deg2rad(50))
        test_getspec_source(meta, source)
    end

    @testset "multi-component sources" begin
        # this test is pretty gnarly
        components = [PointSource("1", Direction(dir"AZEL", 10degrees, 30degrees),
                                  PowerLaw(rand(StokesVector), 10e6, [-rand()])),
                      GaussianSource("2", Direction(dir"AZEL", 10degrees, 30degrees),
                                     PowerLaw(rand(StokesVector), 10e6, [-rand()]),
                                     deg2rad(500/3600), deg2rad(300/3600), deg2rad(50))]
        source = MultiSource("multi", components)
        test_getspec_multi(meta, source)
    end
end

