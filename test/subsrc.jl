let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant

    name,ms = createms(Nant,Nfreq)
    # note that the ms is created with sources already
    # added to the measurement set
    sources = readsources("sources.json")
    subsrc!(ms,sources,TTCal.ConstantBeam())

    @test TTCal.get_corrected_data(ms) ≈ zeros(Complex64,4,Nfreq,Nbase)
end

