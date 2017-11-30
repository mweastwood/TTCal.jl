function old_genvis_equivalent(meta, beam, sky)
    model = zeros(TTCal.JonesMatrix, Nbase(meta), Nfreq(meta))
    frame = ReferenceFrame(meta)
    phase_center = measure(frame, meta.phase_centers[1], dir"ITRF")
    for source in sky.sources
        azel = measure(frame, source.shapes[1].direction, dir"AZEL")
        az = longitude(azel)
        el =  latitude(azel)
        dir = measure(frame, source.shapes[1].direction, dir"ITRF")
        l = dir.x - phase_center.x
        m = dir.y - phase_center.y
        n = dir.z - phase_center.z
        for β = 1:TTCal.Nfreq(meta)
            ν = meta.frequencies[β]
            λ = u"c" / ν
            flux = source.shapes[1].spectrum(ν) |> TTCal.HermitianJonesMatrix
            jones = beam(ν, az, el)
            flux  = TTCal.congruence_transform(jones, flux)
            α = 1
            for antenna1 = 1:Nant(meta), antenna2 = antenna1:Nant(meta)
                r1 = meta.positions[antenna1]
                r2 = meta.positions[antenna2]
                u = (r1.x - r2.x) * u"m" / λ
                v = (r1.y - r2.y) * u"m" / λ
                w = (r1.z - r2.z) * u"m" / λ
                model[α, β] += exp(2π*1im * (u*l + v*m + w*n)) * flux
                α += 1
            end
        end
    end
    model
end

function check_new_and_old_genvis(new, old)
    for β = 1:Nfreq(new)
        α = 1
        for ant1 = 1:Nant(new), ant2 = ant1:Nant(new)
            if !(new[β, 1][ant1, ant2] ≈ old[α, β])
                return false
            end
            α += 1
        end
    end
    true
end

function check_all_unity(dataset)
    for β = 1:Nfreq(dataset)
        α = 1
        for ant1 = 1:Nant(dataset), ant2 = ant1:Nant(dataset)
            if !(dataset[β, 1][ant1, ant2] ≈ one(TTCal.JonesMatrix))
                return false
            end
            α += 1
        end
    end
    true
end

function check_visibilities_match(dataset1, dataset2)
    for β = 1:Nfreq(dataset1)
        α = 1
        for ant1 = 1:Nant(dataset1), ant2 = ant1:Nant(dataset1)
            if !(dataset1[β, 1][ant1, ant2] ≈ dataset2[β, 1][ant1, ant2])
                return false
            end
            α += 1
        end
    end
    true
end

@testset "genvis.jl" begin
    Nant  = 5
    Nfreq = 10
    name, ms = createms(Nant, Nfreq)
    metadata = TTCal.Metadata(ms)
    beam = TTCal.ConstantBeam()
    sky  = readsky("sources.json")

    #@testset "geometric delays" begin
    #    # check that the near field and far field routines give equal geometric delays
    #    # in the far field
    #    position = TTCal.position(meta)
    #    vector = [position.x; position.y; position.z]
    #    farfield_vector  = vector / norm(vector)
    #    nearfield_vector = vector * 1e5
    #    farfield_direction = Direction(dir"ITRF",  farfield_vector[1],  farfield_vector[2],  farfield_vector[3])
    #    nearfield_position =  Position(pos"ITRF", nearfield_vector[1], nearfield_vector[2], nearfield_vector[3])
    #    farfield_delays  = TTCal.geometric_delays(meta.antennas, farfield_direction, meta.phase_center)
    #    nearfield_delays = TTCal.geometric_delays(meta.antennas, nearfield_position, meta.phase_center)
    #    @test vecnorm(farfield_delays - nearfield_delays)/vecnorm(farfield_delays) < 1e-8
    #end

    @testset "point sources" begin
        # test that genvis agrees with the above reference implementation
        new = genvis(metadata, beam, sky)
        old = old_genvis_equivalent(metadata, beam, sky)
        @test check_new_and_old_genvis(new, old)

        # a source at the phase center should have unity visibilities
        sources = [TTCal.Source("Cristiano Ronaldo",
                                TTCal.Point(metadata.phase_centers[1],
                                            TTCal.PowerLaw(1, 0, 0, 0, 10*u"MHz", [0.0])))]
        dataset = genvis(metadata, beam, TTCal.SkyModel(sources))
        @test check_all_unity(dataset)
    end

    @testset "extended sources" begin
        # the general strategy here will be to check that the visibilities match those
        # of a point source in the appropriate limits
        point = TTCal.Source("point",
                             TTCal.Point(Direction(dir"AZEL", 45*u"°", 45*u"°"),
                                         TTCal.PowerLaw(4, 3, 2, 1, 10*u"MHz", [1.0])))
        gaussian = TTCal.Source("gaussian",
                                TTCal.Gaussian(Direction(dir"AZEL", 45*u"°", 45*u"°"),
                                               TTCal.PowerLaw(4, 3, 2, 1, 10*u"MHz", [1.0]),
                                               0., 0., 0.))
        #disk = DiskSource("disk", Direction(dir"AZEL", 45degrees, 45degrees),
        #                          PowerLaw(4, 3, 2, 1, 10e6, [1.0]), 0)
        #shapelet = ShapeletSource("shapelet", Direction(dir"AZEL", 45degrees, 45degrees),
        #                                      PowerLaw(4, 3, 2, 1, 10e6, [1.0]), 0, [1.0])
        point_visibilities = genvis(metadata, beam, TTCal.SkyModel(point))
        gaussian_visibilities = genvis(metadata, beam, TTCal.SkyModel(gaussian))
        #disk_visibilities = genvis(metadata, disk)
        #shapelet_visibilities = genvis(metadata, shapelet)
        @test check_visibilities_match(point_visibilities, gaussian_visibilities)
        #threshold = eps(Float64)*vecnorm(point_visibilities.data)
        #@test vecnorm(point_visibilities.data - gaussian_visibilities.data) < threshold
        #@test vecnorm(point_visibilities.data - disk_visibilities.data) < threshold
        #@test vecnorm(point_visibilities.data - shapelet_visibilities.data) < threshold

        # the first shapelet component is a Gaussian so we can check that the Gaussian and
        # shapelet visibilities match in this limit
        #let
        #    β = deg2rad(1) # the shapelet scale parameter
        #    fwhm = asin(2sqrt(2log(2))*β) # the Gaussian full-width-half-maximum
        #    # note that the `asin` comes about because of a small angle approximation that is
        #    # made in the shapelet `genvis` but not Gaussian `genvis`
        #    gaussian = GaussianSource("gaussian", Direction(dir"AZEL", 45degrees, 45degrees),
        #                                          PowerLaw(4, 3, 2, 1, 10e6, [1.0]), fwhm, fwhm, 0)
        #    shapelet = ShapeletSource("shapelet", Direction(dir"AZEL", 45degrees, 45degrees),
        #                                          PowerLaw(4, 3, 2, 1, 10e6, [1.0]), β, [1.0])
        #    gaussian_visibilities = genvis(meta, gaussian)
        #    shapelet_visibilities = genvis(meta, shapelet)
        #    threshold = eps(Float64)*vecnorm(gaussian_visibilities.data)
        #    @test vecnorm(shapelet_visibilities.data - gaussian_visibilities.data) < threshold
        #end
    end

    #@testset "multi component sources" begin
    #    # a multi-component source should get the same visibilities as each of those sources together
    #    sources = [PointSource("point", Direction(dir"AZEL", 45degrees, 45degrees),
    #                                    PowerLaw(4, 3, 2, 1, 10e6, [1.0])),
    #               GaussianSource("gaussian", Direction(dir"AZEL", 45degrees, 45degrees),
    #                                          PowerLaw(4, 3, 2, 1, 10e6, [1.0]),
    #                                          deg2rad(500/3600), deg2rad(300/3600), deg2rad(50))]
    #    multi = MultiSource("multi", sources)
    #    vis1 = genvis(meta, sources)
    #    vis2 = genvis(meta, multi)
    #    @test vecnorm(vis1.data - vis2.data) < eps(Float64) * vecnorm(vis1.data)
    #end

    Tables.delete(ms)
end

