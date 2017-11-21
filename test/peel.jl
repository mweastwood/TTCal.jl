@testset "peel.jl" begin
    Nant = 10
    Nfreq = 2

    name,ms = createms(Nant, Nfreq)
    meta = Metadata(ms)
    beam = ConstantBeam()
    sources = readsources("sources.json")

    @testset "peeling" begin
        let
            visibilities = genvis(meta, sources)
            δ = sqrt(eps(Float64))*vecnorm(visibilities.data)
            solutions = peel!(visibilities, meta, beam, sources, maxiter=100, tolerance=eps(Float64))
            for solution in solutions
                @test vecnorm(solution.jones - ones(DiagonalJonesMatrix,Nant,Nfreq)) < 1e-6
            end
            @test vecnorm(visibilities.data) < δ
            @test !any(visibilities.flags)
        end

        let
            visibilities = genvis(meta, sources)
            δ = sqrt(eps(Float64))*vecnorm(visibilities.data)
            TTCal.write(ms, "DATA", visibilities)
            Tables.unlock(ms)
            TTCal.main(["peel", name, "sources.json", "--minuvw", "0",
                        "--peeliter", "3", "--maxiter", "100", "--tolerance","$(eps(Float64))",
                        "--beam", "constant"])
            Tables.lock(ms)
            visibilities = TTCal.read(ms, "DATA")
            @test vecnorm(visibilities.data) < 10δ
            @test !any(visibilities.flags)
        end
    end

    @testset "shaving" begin
        let
            visibilities = genvis(meta, sources)
            δ = sqrt(eps(Float64))*vecnorm(visibilities.data)
            solutions = shave!(visibilities, meta, beam, sources, maxiter=100, tolerance=eps(Float64))
            for solution in solutions
                @test vecnorm(solution.jones - ones(DiagonalJonesMatrix,Nant,1)) < 1e-6
            end
            @test vecnorm(visibilities.data) < δ
            @test !any(visibilities.flags)
        end

        let
            visibilities = genvis(meta, sources)
            δ = sqrt(eps(Float64))*vecnorm(visibilities.data)
            TTCal.write(ms, "DATA", visibilities)
            Tables.unlock(ms)
            TTCal.main(["shave", name, "sources.json", "--minuvw", "0",
                        "--peeliter", "3", "--maxiter", "100", "--tolerance","$(eps(Float64))",
                        "--beam", "constant"])
            Tables.lock(ms)
            visibilities = TTCal.read(ms, "DATA")
            @test vecnorm(visibilities.data) < 10δ
            @test !any(visibilities.flags)
        end
    end
end

