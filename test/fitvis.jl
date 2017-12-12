@testset "fitvis.jl" begin
    positions = [Position(pos"ITRF", 1e2randn(), 1e2randn(), 1e2randn()) for n = 1:10]
    metadata = TTCal.Metadata([50.0u"MHz", 74.3u"MHz", 93.1u"MHz"],
                              [Epoch(epoch"UTC", (2017-1858)*u"yr")],
                              positions,
                              [Direction(dir"AZEL", 0u"°", 90u"°")])
    tolerance = eps(Float64)
    spectrum  = TTCal.PowerLaw(1, 0, 0, 0, 100u"MHz", [-0.5])
    direction = Direction(dir"AZEL", 10u"°", 30u"°")
    itrf = measure(ReferenceFrame(metadata), direction, dir"ITRF")
    ra   = longitude(itrf)
    dec  =  latitude(itrf)
    δra  = uconvert(u"rad", 10randn()*u"arcminute")
    δdec = uconvert(u"rad", 10randn()*u"arcminute")
    itrf′ = Direction(dir"ITRF", ra+δra, dec+δdec)

    @testset "point sources" begin
        before = TTCal.Source("A", TTCal.Point(itrf,  spectrum))
        after  = TTCal.Source("B", TTCal.Point(itrf′, spectrum))
        data = genvis(metadata, TTCal.ConstantBeam(), after, polarization=TTCal.Dual)
        fit  = fitvis(data, before, tolerance=tolerance)
        @test fit.shapes[1].direction ≈ after.shapes[1].direction
    end

    @testset "gaussian sources" begin
        a = deg2rad(500/3600)
        b = deg2rad(300/3600)
        c = deg2rad(50)
        before = TTCal.Source("A", TTCal.Gaussian(itrf,  spectrum, a, b, c))
        after  = TTCal.Source("B", TTCal.Gaussian(itrf′, spectrum, a, b, c))
        data = genvis(metadata, TTCal.ConstantBeam(), after, polarization=TTCal.Dual)
        fit  = fitvis(data, before, tolerance=tolerance)
        @test fit.shapes[1].direction ≈ after.shapes[1].direction
    end

#    @testset "multi-component sources" begin
#        components = [PointSource("1", Direction(dir"AZEL", 10degrees, 30degrees),
#                                  PowerLaw(1, 0, 0, 0, 10e6, [-rand()])),
#                      GaussianSource("2", Direction(dir"AZEL", 10degrees, 30degrees),
#                                     PowerLaw(1, 0, 0, 0, 10e6, [-rand()]),
#                                     deg2rad(500/3600), deg2rad(300/3600), deg2rad(50))]
#        source = MultiSource("multi", components)
#        direction = measure(frame, components[1].direction, dir"J2000")
#        data = genvis(meta, source)
#        ra   = longitude(direction)
#        dec  =  latitude(direction)
#        ra  += 0.1 * π/180 * randn()
#        dec += 0.1 * π/180 * randn()
#        for component in components
#            component.direction = Direction(dir"J2000", ra*radians, dec*radians)
#        end
#        fitdirection = fitvis(data, meta, source, tolerance = tolerance)
#        @test fitdirection ≈ direction
#    end
end

