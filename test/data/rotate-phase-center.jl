@testset "rotate-phase-center.jl" begin
    metadata = TTCal.Metadata([50.0u"MHz", 74.3u"MHz", 93.1u"MHz"],
                              [Epoch(epoch"UTC", (2017-1858)*u"yr")],
                              [Position(pos"ITRF", 100randn(), 100randn(), 100randn()) for n = 1:3],
                              [Direction(dir"AZEL", 0u"°", 90u"°")])

    source = TTCal.Source("A thing",
                          TTCal.Point(Direction(dir"J2000", "0h12m34s", "+12d34m56s"),
                                      TTCal.PowerLaw(1, 0, 0, 0, 100u"MHz", [0.0])))
    dataset = genvis(metadata, TTCal.ConstantBeam(), source)
    rotate_phase_center!(dataset, Direction(dir"J2000", "0h12m34s", "+12d34m56s"))

    for frequency = 1:Nfreq(dataset)
        visibilities = dataset[frequency, 1]
        for ant1 = 1:Nant(dataset), ant2 = ant1:Nant(dataset)
            @test visibilities[ant1, ant2] ≈ one(TTCal.JonesMatrix)
        end
    end
end

