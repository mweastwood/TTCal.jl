let I0 = 100, Q0 = -50, U0 = 40, V0 = -10, ν0 = 123, index = [1.0]
    spec = Spectrum(I0,Q0,U0,V0,ν0,index)

    frequency = 123.0
    @test norm(spec(frequency) - StokesVector(I0,Q0,U0,V0)) < eps(Float64)
    frequency = 1230.0
    @test norm(spec(frequency) - 10*StokesVector(I0,Q0,U0,V0)) < eps(Float64)
end

let
    frame   = ReferenceFrame()
    sources = [Source("S1",Point("P1",Direction(dir"AZEL","0.0deg","+45.0deg"),
                           Spectrum(1.0,2.0,3.0,4.0,10e6,[0.0]))),
               Source("S2",Point("P2",Direction(dir"AZEL","0.0deg","-45.0deg"),
                           Spectrum(1.0,2.0,3.0,4.0,10e6,[0.0])))]
    @test TTCal.isabovehorizon(frame,sources[1]) == true
    @test TTCal.isabovehorizon(frame,sources[2]) == false

    sources′ = TTCal.abovehorizon(frame,sources)
    @test length(sources′) == 1
    @test sources′[1].name == "S1"
end

let
    sources = readsources("sources.json")

    cas_a = sources[1]
    point = cas_a.components[1]
    @test cas_a.name == "Cas A"
    @test point.name == "Cas A"
    @test point.direction == Direction(dir"J2000","23h23m24s","58d48m54s")
    @test point.spectrum.stokes.I == 555904.26
    @test point.spectrum.stokes.Q == 0
    @test point.spectrum.stokes.U == 0
    @test point.spectrum.stokes.V == 0
    @test point.spectrum.ν0 == 1e6
    @test point.spectrum.spectral_index == [-0.770]

    cyg_a = sources[2]
    point = cyg_a.components[1]
    @test cyg_a.name == "Cyg A"
    @test point.name == "Cyg A"
    @test point.direction == Direction(dir"J2000","19h59m28.35663s","+40d44m02.0970s")
    @test point.spectrum.stokes.I == 49545.02
    @test point.spectrum.stokes.Q == 0
    @test point.spectrum.stokes.U == 0
    @test point.spectrum.stokes.V == 0
    @test point.spectrum.ν0 == 1e6
    @test point.spectrum.spectral_index == [+0.085,-0.178]
end

let
    sources1 = [Source("S1",Point("P1",Direction(dir"J2000","1h","0d"),
                            Spectrum(1.0,2.0,3.0,4.0,10e6,[0.0]))),
                Source("S2",Point("P2",Direction(dir"J2000","2h","0d"),
                            Spectrum(1.0,2.0,3.0,4.0,10e6,[0.0]))),
                Source("S3",Point("P3",Direction(dir"J2000","3h","0d"),
                            Spectrum(1.0,2.0,3.0,4.0,10e6,[0.0])))]

    filename = tempname()*".json"
    writesources(filename,sources1)
    sources2 = readsources(filename)
    for i in eachindex(sources1,sources2)
        source1 = sources1[i]
        source2 = sources2[i]
        point1  = source1.components[1]
        point2  = source2.components[1]
        @test source1.name == source2.name
        @test point1.name == point2.name
        @test point1.direction == point2.direction
        @test point1.spectrum.stokes.I == point2.spectrum.stokes.I
        @test point1.spectrum.stokes.Q == point2.spectrum.stokes.Q
        @test point1.spectrum.stokes.U == point2.spectrum.stokes.U
        @test point1.spectrum.stokes.V == point2.spectrum.stokes.V
        @test point1.spectrum.ν0 == point2.spectrum.ν0
        @test point1.spectrum.spectral_index == point2.spectrum.spectral_index
    end
end

let
    frame = ReferenceFrame()
    t = (2015.-1858.)*365. # a rough, current Julian date
    set!(frame,Epoch(epoch"UTC",t * days))
    set!(frame,observatory("OVRO_MMA"))
    phase_dir = Direction(dir"J2000","19h59m28.35663s","+40d44m02.0970s")

    for iteration = 1:5
        l = 2rand()-1
        m = 2rand()-1
        while hypot(l,m) > 1
            l = 2rand()-1
            m = 2rand()-1
        end
        dir = TTCal.undo_direction_cosines(phase_dir,l,m)
        l′,m′ = TTCal.direction_cosines(phase_dir,dir)
        @test l ≈ l′
        @test m ≈ m′
    end

    # z should be negative in undo_direction_cosines
    l,m = (-0.26521340920368575,-0.760596242177856)
    dir = TTCal.undo_direction_cosines(phase_dir,l,m)
    l′,m′ = TTCal.direction_cosines(phase_dir,dir)
    @test l ≈ l′
    @test m ≈ m′
end

