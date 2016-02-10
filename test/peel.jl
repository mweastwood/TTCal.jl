@testset "peel.jl" begin
    Nant = 10
    Nfreq = 2

    name,ms = createms(Nant, Nfreq)
    meta = collect_metadata(ms, ConstantBeam())
    sources = readsources("sources.json")
    visibilities = genvis(meta, sources)

    solutions = peel!(GainCalibration, visibilities, meta, sources, maxiter=100, tolerance=eps(Float64))
    for solution in solutions
        @test vecnorm(solution.jones - ones(DiagonalJonesMatrix,Nant,Nfreq)) < 1e-6
    end
end

