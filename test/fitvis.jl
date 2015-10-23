let
    Nant = 10
    Nfreq = 2
    name,ms = createms(Nant,Nfreq)

    # Regenerate the data column with only one source
    # to prevent sidelobe contamination in the tests
    sources = readsources("sources.json")[1:1]
    ms.table["DATA"] = genvis(ms.frame,ms.phase_direction,sources,
                              TTCal.ConstantBeam(),ms.u,ms.v,ms.w,ms.ν)

    # Perturb the source position
    az,el = TTCal.azel(ms.frame,sources[1])
    lold,mold = TTCal.lm(ms.phase_direction,sources[1])
    az += 0.1 * π/180
    el += 0.1 * π/180
    sources[1].dir = Direction(dir"AZEL",Quantity(az,"rad"),
                                         Quantity(el,"rad"))

    l,m = fitvis(ms,sources,maxiter = 100, tolerance = sqrt(eps(Float64)))

    @test l[1] ≈ lold
    @test m[1] ≈ mold
end

