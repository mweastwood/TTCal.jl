@testset "peel.jl" begin
    Nant = 10
    Nfreq = 2

    name,ms = createms(Nant, Nfreq)
    meta = collect_metadata(ms, ConstantBeam())
    sources = readsources("sources.json")

    @testset "peeling" begin
        let
            visibilities = genvis(meta, sources)
            δ = sqrt(eps(Float64))*vecnorm(visibilities.data)
            solutions = peel!(visibilities, meta, sources, maxiter=100, tolerance=eps(Float64))
            for solution in solutions
                @test vecnorm(solution.jones - ones(DiagonalJonesMatrix,Nant,Nfreq)) < 1e-6
            end
            @test vecnorm(visibilities.data) < δ
            @test !any(visibilities.flags)
        end

        let
            visibilities = genvis(meta, sources)
            δ = sqrt(eps(Float64))*vecnorm(visibilities.data)
            TTCal.set_data!(ms, visibilities)
            unlock(ms)
            TTCal.main(["peel","--input",name,"--sources","sources.json","--minuvw","0",
                        "--peeliter","3","--maxiter","100","--tolerance","$(eps(Float64))","--beam","constant"])
            lock(ms)
            visibilities = TTCal.get_corrected_data(ms)
            @test vecnorm(visibilities.data) < 10δ
            @test !any(visibilities.flags)
        end
    end

    @testset "shaving" begin
        let
            visibilities = genvis(meta, sources)
            δ = sqrt(eps(Float64))*vecnorm(visibilities.data)
            solutions = shave!(visibilities, meta, sources, maxiter=100, tolerance=eps(Float64))
            for solution in solutions
                @test vecnorm(solution.jones - ones(DiagonalJonesMatrix,Nant,1)) < 1e-6
            end
            @test vecnorm(visibilities.data) < δ
            @test !any(visibilities.flags)
        end

        let
            visibilities = genvis(meta, sources)
            δ = sqrt(eps(Float64))*vecnorm(visibilities.data)
            TTCal.set_data!(ms, visibilities)
            unlock(ms)
            TTCal.main(["shave","--input",name,"--sources","sources.json","--minuvw","0",
                        "--peeliter","3","--maxiter","100","--tolerance","$(eps(Float64))","--beam","constant"])
            lock(ms)
            visibilities = TTCal.get_corrected_data(ms)
            @test vecnorm(visibilities.data) < 10δ
            @test !any(visibilities.flags)
        end
    end
end

