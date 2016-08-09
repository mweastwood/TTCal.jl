@testset "beammodel.jl" begin
    @testset "constant beam" begin
        beam = ConstantBeam()
        for i = 1:3
            ν  = rand() * 100e6
            az = rand() * 2π
            el = rand() * π/2
            @test beam(ν,az,el) == one(JonesMatrix)
        end
    end

    @testset "sine beam" begin
        α = 2.0
        beam = SineBeam(α)
        for i = 1:3
            ν  = rand() * 100e6
            az = rand() * 2π
            el = rand() * π/2
            M = beam(ν,az,el) |> MuellerMatrix
            @test Matrix(M) ≈ sin(el)^α * eye(4)
        end

        beam1 = SineBeam()
        beam2 = SineBeam(1.6)
        for i = 1:3
            ν  = rand() * 100e6
            az = rand() * 2π
            el = rand() * π/2
            @test beam1(ν,az,el) == beam2(ν,az,el)
        end

    end

    @testset "memo 178 beam" begin
        beam = Memo178Beam()
        @test beam(45e6,0,π/2) == one(JonesMatrix)
    end

    @testset "Zernike beam" begin
        # These coefficients are from a fit to the OVRO LWA Stokes I beam at 67.920 MHz.
        # One of the constraints on the fit was that the beam is unity at zenith.
        coeff = [ 0.5438736462155408,
                 -0.451837602076611,
                 -0.01674201596093591,
                 -0.0403607714687312,
                 -0.032646170124955354,
                  0.04700116243381884,
                 -0.011618155119718811,
                 -0.0008121622545099635,
                 -0.016239776458086115]
        beam = ZernikeBeam(coeff)
        @test norm(beam(67.920e6,0,π/2) - one(JonesMatrix)) < 1e-5
    end
end

