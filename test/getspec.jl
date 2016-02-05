@testset "getspec.jl" begin
    Nant = 10
    Nfreq = 2
    name,ms = createms(Nant,Nfreq)
    meta = TTCal.collect_metadata(ms, ConstantBeam())

    # Generate the data column with only one source
    # to prevent sidelobe contamination in the tests
    sources = readsources("sources.json")[1:1]
    data = genvis(meta, sources)

    dir = sources[1].direction
    stokes_flux = sources[1].spectrum(meta.channels)
    measured_flux = getspec(data, meta, dir)
    for β = 1:Nfreq
        linear_flux = TTCal.linear(stokes_flux[β])
        @test linear_flux.xx ≈ measured_flux[β].xx
        @test linear_flux.xy ≈ measured_flux[β].xy
        @test linear_flux.yy ≈ measured_flux[β].yy
    end

    # try again with a polarized source
    source = PointSource("FRB", meta.phase_center, PowerLaw(rand(StokesVector),10e6,[1.0]))
    data = genvis(meta, [source])

    dir = source.direction
    stokes_flux = source.spectrum(meta.channels)
    measured_flux = getspec(data, meta, dir)
    for β = 1:Nfreq
        linear_flux = TTCal.linear(stokes_flux[β])
        @test linear_flux.xx ≈ measured_flux[β].xx
        @test linear_flux.xy ≈ measured_flux[β].xy
        @test linear_flux.yy ≈ measured_flux[β].yy
    end
end

