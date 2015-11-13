let
    Nant = 10
    Nfreq = 2
    name,ms = createms(Nant,Nfreq)

    # Generate the data column with only one source
    # to prevent sidelobe contamination in the tests
    source = readsources("sources.json")[1]
    point  = source.components[1]
    ms.table["DATA"] = genvis(ms, source, TTCal.ConstantBeam())

    stokes_flux = point.spectrum(ms.ν)
    measured_flux = getspec(ms,point.direction)
    for β = 1:Nfreq
        linear_flux = TTCal.linear(stokes_flux[β])
        @test linear_flux.xx ≈ measured_flux[β].xx
        @test linear_flux.xy ≈ measured_flux[β].xy
        @test linear_flux.yy ≈ measured_flux[β].yy
    end

    # try again with a polarized source
    source = Source("FRB",
                    Point("from z=1",
                          ms.phase_direction,
                          Spectrum(rand(StokesVector),10e6,[1.0])))
    point  = source.components[1]
    ms.table["DATA"] = genvis(ms, source, TTCal.ConstantBeam())

    stokes_flux = point.spectrum(ms.ν)
    measured_flux = getspec(ms,point.direction)
    for β = 1:Nfreq
        linear_flux = TTCal.linear(stokes_flux[β])
        @test isapprox(linear_flux.xx, measured_flux[β].xx,atol=1e-6)
        @test isapprox(linear_flux.xy, measured_flux[β].xy,atol=1e-6)
        @test isapprox(linear_flux.yy, measured_flux[β].yy,atol=1e-6)
    end
end

