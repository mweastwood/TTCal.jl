let
    ϕ = linspace(1,10,100)
    fringe_naive = exp(1im*ϕ)
    fringe = TTCal.fringepattern(ϕ[1],ϕ[2]-ϕ[1],length(ϕ))
    @test isapprox(fringe,fringe_naive,rtol=10eps(Float32))
end

