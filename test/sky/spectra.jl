@testset "spectra.jl" begin
    @testset "power law" begin
        I0 = 100; Q0 = -50; U0 = 40; V0 = -10; ν0 = 123u"Hz"; index = [1.0]

        spec = TTCal.PowerLaw(I0, Q0, U0, V0, ν0, index)
        @test spec( 123.0*u"Hz") ≈    TTCal.StokesVector(I0, Q0, U0, V0)
        @test spec(1230.0*u"Hz") ≈ 10*TTCal.StokesVector(I0, Q0, U0, V0)

        spec = TTCal.PowerLaw(-I0, -Q0, -U0, -V0, ν0, index)
        @test spec( 123.0*u"Hz") ≈    TTCal.StokesVector(-I0, -Q0, -U0, -V0)
        @test spec(1230.0*u"Hz") ≈ 10*TTCal.StokesVector(-I0, -Q0, -U0, -V0)
    end

    #e = one(TTCal.StokesVector)
    #spec = RFISpectrum([1e6, 2e6, 3e6], [1e, 2e, 3e])
    #@test spec(1e6) == 1one(TTCal.StokesVector)
    #@test spec(2e6) == 2one(TTCal.StokesVector)
    #@test spec(3e6) == 3one(TTCal.StokesVector)
end

