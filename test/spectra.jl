@testset "spectra.jl" begin
    I0 = 100; Q0 = -50; U0 = 40; V0 = -10; ν0 = 123; index = [1.0]
    spec = PowerLaw(I0,Q0,U0,V0,ν0,index)
    @test norm(spec(123.0) - StokesVector(I0,Q0,U0,V0)) < eps(Float64)
    @test norm(spec(1230.0) - 10*StokesVector(I0,Q0,U0,V0)) < eps(Float64)

    e = one(StokesVector)
    spec = RFISpectrum([1e6, 2e6, 3e6], [1e, 2e, 3e])
    @test spec(1e6) == 1one(StokesVector)
    @test spec(2e6) == 2one(StokesVector)
    @test spec(3e6) == 3one(StokesVector)
end

