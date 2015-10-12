let
    Nant = 10
    Nfreq = 2
    name,ms = createms(Nant,Nfreq)

    # Regenerate the data column with only one source
    # to prevent sidelobe contamination in the tests
    sources = readsources("sources.json")[1:1]
    ms.table["DATA"] = genvis(ms.frame,ms.phase_direction,sources,
                              TTCal.ConstantBeam(),ms.u,ms.v,ms.w,ms.ν)

    flux = TTCal.flux(sources[1],ms.ν)
    xx,xy,yx,yy = getspec(ms,TTCal.direction(sources[1]))
    @test xx ≈ flux
    @test yy ≈ flux
end

