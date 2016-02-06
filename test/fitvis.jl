@testset "fitvis.jl" begin
    Nant = 10
    Nfreq = 2
    name,ms = createms(Nant,Nfreq)
    meta = TTCal.collect_metadata(ms, ConstantBeam())
    sources = readsources("sources.json")[1:1]
    data = genvis(meta, sources)
    # Perturb the source position
    direction = sources[1].direction
    ra   = longitude(direction)
    dec  =  latitude(direction)
    ra  += 0.1 * π/180 * randn()
    dec += 0.1 * π/180 * randn()
    newdirection = Direction(dir"J2000", ra*radians, dec*radians)
    fitdirection = fitvis(data, meta, newdirection, maxiter = 100, tolerance = sqrt(eps(Float64)))
    @test fitdirection ≈ direction
end

