let
    Nant = 10
    Nfreq = 2
    name,ms = createms(Nant,Nfreq)

    # Regenerate the data column with only one source
    # to prevent sidelobe contamination in the tests
    sources = readsources("sources.json")[1:1]
    point   = sources[1].components[1]
    ms.table["DATA"] = genvis(ms,sources,TTCal.ConstantBeam())

    # Perturb the source position
    az,el = TTCal.local_azel(ms.frame,point)
    lold,mold = TTCal.direction_cosines(ms.phase_direction,point.direction)
    az += 0.1 * π/180
    el += 0.1 * π/180
    dir = Direction(dir"AZEL",Quantity(az,"rad"),Quantity(el,"rad"))

    l,m = fitvis(ms,dir,maxiter = 100, tolerance = sqrt(eps(Float64)))

    @test l[1] ≈ lold
    @test m[1] ≈ mold
end

