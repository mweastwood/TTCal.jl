import TTCal: name, direction

let
    sources = readsources("sources.json")

    cas_a = sources[1]
    @test name(cas_a) == "Cas A"
    dir1 = direction(cas_a)
    dir2 = Direction(Measures.J2000,ra"23h23m24s",dec"58d48m54s")
    @test Measures.reference(dir1) == Measures.reference(dir2)
    @test longitude(dir1) ≈ longitude(dir2)
    @test latitude(dir1) ≈ latitude(dir2)
    @test cas_a.I == 555904.26
    @test cas_a.Q == 0.0
    @test cas_a.U == 0.0
    @test cas_a.V == 0.0
    @test cas_a.reffreq == 1.0e6
    @test cas_a.index == [-0.770]

    cyg_a = sources[2]
    @test name(cyg_a) == "Cyg A"
    dir1 = direction(cyg_a)
    dir2 = Direction(Measures.J2000,ra"19h59m28.35663s",dec"+40d44m02.0970s")
    @test Measures.reference(dir1) == Measures.reference(dir2)
    @test longitude(dir1) ≈ longitude(dir2)
    @test latitude(dir1) ≈ latitude(dir2)
    @test cyg_a.I == 49545.02
    @test cyg_a.Q == 0.0
    @test cyg_a.U == 0.0
    @test cyg_a.V == 0.0
    @test cyg_a.reffreq == 1.0e6
    @test cyg_a.index == [+0.085,-0.178]
end

let
    sources1 = [PointSource("S1",Direction(Measures.J2000,ra"1h",dec"0d"),
                            1.0,2.0,3.0,4.0,10e6,[0.0]),
                PointSource("S2",Direction(Measures.J2000,ra"2h",dec"0d"),
                            1.0,2.0,3.0,4.0,10e6,[0.0]),
                PointSource("S3",Direction(Measures.J2000,ra"3h",dec"0d"),
                            1.0,2.0,3.0,4.0,10e6,[0.0])]
    filename = tempname()
    writesources(filename,sources1)
    sources2 = readsources(filename)
    for i in eachindex(sources1,sources2)
        @test name(sources1[i]) == name(sources2[i])
        dir1 = direction(sources1[i])
        dir2 = direction(sources2[i])
        @test Measures.reference(dir1) == Measures.reference(dir2)
        @test longitude(dir1) ≈ longitude(dir2)
        @test latitude(dir1) ≈ latitude(dir2)
        @test sources1[i].I == sources2[i].I
        @test sources1[i].Q == sources2[i].Q
        @test sources1[i].U == sources2[i].U
        @test sources1[i].V == sources2[i].V
        @test sources1[i].reffreq == sources2[i].reffreq
        @test sources1[i].index == sources2[i].index
    end
end

let
    source = PointSource("S",Direction(Measures.J2000,ra"1h",dec"0d"),
                         1.0,2.0,3.0,4.0,10e6,[0.0])
    @test name(source) == "S"
    @test repr(source) == "S"
    @test TTCal.flux(source,10e6) == 1.0
    @test TTCal.flux(source,20e6) == 1.0
end

let
    frame = ReferenceFrame()
    t = (2015.-1858.)*365.*24.*60.*60. # a rough, current Julian date (in seconds)
    set!(frame,Epoch(Measures.UTC,Quantity(t,Second)))
    set!(frame,observatory("OVRO_MMA"))
    phase_dir = Direction(Measures.J2000,ra"19h59m28.35663s",dec"+40d44m02.0970s")

    for iteration = 1:5
        l = 2rand()-1
        m = 2rand()-1
        while hypot(l,m) > 1
            l = 2rand()-1
            m = 2rand()-1
        end
        dir = TTCal.lm2dir(phase_dir,l,m)
        l′,m′ = TTCal.dir2lm(frame,phase_dir,dir)
        @test l ≈ l′
        @test m ≈ m′
    end

    # z should be negative in lm2dir
    l,m = (-0.26521340920368575,-0.760596242177856)
    dir = TTCal.lm2dir(phase_dir,l,m)
    l′,m′ = TTCal.dir2lm(frame,phase_dir,dir)
    @test l ≈ l′
    @test m ≈ m′
end

