@testset "fitvis.jl" begin
    Nant = 20
    Nfreq = 2
    name, ms = createms(Nant, Nfreq)
    meta = Metadata(ms)
    beam = ConstantBeam()
    frame = TTCal.reference_frame(meta)
    tolerance = eps(Float64)/1e5

    @testset "point sources" begin
        sources = readsources("sources.json")[1:1]
        data = genvis(meta, sources)
        direction = sources[1].direction
        ra   = longitude(direction)
        dec  =  latitude(direction)
        ra  += 0.1 * π/180 * randn()
        dec += 0.1 * π/180 * randn()
        newdirection = Direction(dir"J2000", ra*radians, dec*radians)
        fitdirection = fitvis(data, meta, newdirection, tolerance = tolerance)
        @test fitdirection ≈ direction
    end

    @testset "gaussian sources" begin
        spectrum = PowerLaw(rand(StokesVector), 10e6, [-rand()])
        direction = measure(frame, Direction(dir"AZEL", 10degrees, 30degrees), dir"J2000")
        source = GaussianSource("test", direction, spectrum,
                                        deg2rad(500/3600), deg2rad(300/3600), deg2rad(50))

        data = genvis(meta, source)
        ra   = longitude(direction)
        dec  =  latitude(direction)
        ra  += 0.1 * π/180 * randn()
        dec += 0.1 * π/180 * randn()
        source.direction = Direction(dir"J2000", ra*radians, dec*radians)
        fitdirection = fitvis(data, meta, source, tolerance = tolerance)
        @test fitdirection ≈ direction
    end

    @testset "multi-component sources" begin
        components = [PointSource("1", Direction(dir"AZEL", 10degrees, 30degrees),
                                  PowerLaw(1, 0, 0, 0, 10e6, [-rand()])),
                      GaussianSource("2", Direction(dir"AZEL", 10degrees, 30degrees),
                                     PowerLaw(1, 0, 0, 0, 10e6, [-rand()]),
                                     deg2rad(500/3600), deg2rad(300/3600), deg2rad(50))]
        source = MultiSource("multi", components)
        direction = measure(frame, components[1].direction, dir"J2000")
        data = genvis(meta, source)
        ra   = longitude(direction)
        dec  =  latitude(direction)
        ra  += 0.1 * π/180 * randn()
        dec += 0.1 * π/180 * randn()
        for component in components
            component.direction = Direction(dir"J2000", ra*radians, dec*radians)
        end
        fitdirection = fitvis(data, meta, source, tolerance = tolerance)
        @test fitdirection ≈ direction
    end
end

