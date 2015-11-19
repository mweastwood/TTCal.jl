let
    Nant = 10
    Nfreq = 2

    name,ms = createms(Nant,Nfreq)
    sources = readsources("sources.json")
    ms.table["DATA"] = genvis(ms,sources,ConstantBeam())

    solutions = peel!(GainCalibration,ms,sources,ConstantBeam(),maxiter=100,tolerance=eps(Float64))
    for solution in solutions
        @test vecnorm(solution.jones - ones(DiagonalJonesMatrix,Nant,Nfreq)) < 1e-6
    end
end

