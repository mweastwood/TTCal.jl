let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant

    name,ms = createms(Nant,Nfreq)
    sources = readsources("sources.json")
    ms.table["DATA"] = genvis(ms,sources,TTCal.ConstantBeam())
    subsrc!(ms,sources,TTCal.ConstantBeam())

    @test TTCal.get_corrected_data(ms) ≈ zeros(Complex64,4,Nfreq,Nbase)
end

let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant

    name,ms = createms(Nant,Nfreq)
    source  = readsources("sources.json")[1]
    ms.table["DATA"] = genvis(ms,source,TTCal.ConstantBeam())
    subsrc!(ms,source.components[1].direction)

    @test isapprox(TTCal.get_corrected_data(ms),zeros(Complex64,4,Nfreq,Nbase),atol=1e-1)
    # note the rough tolerance here has to do with the accuracy of the inverse Fourier
    # transform used by getspec
end

