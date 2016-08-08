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

    @testset "Zernike polynomials" begin
        @test TTCal.zernike(0,  0, 0.5, π/2) == 1
        @test TTCal.zernike(0,  0, 1.0, π/4) == 1

        R11(ρ) = ρ
        @test TTCal.zernike(1, +1, 0.5, π/3) ≈ R11(0.5)*cos(π/3)
        @test TTCal.zernike(1, +1, 1.0, π/6) ≈ R11(1.0)*cos(π/6)
        @test TTCal.zernike(1, -1, 0.5, π/3) ≈ R11(0.5)*sin(π/3)
        @test TTCal.zernike(1, -1, 1.0, π/6) ≈ R11(1.0)*sin(π/6)

        R20(ρ) = 2ρ^2 - 1
        R22(ρ) = ρ^2
        @test TTCal.zernike(2,  0, 0.3, π/3) ≈ R20(0.3)
        @test TTCal.zernike(2, +2, 0.3, π/3) ≈ R22(0.3)*cos(2π/3)
        @test TTCal.zernike(2, -2, 0.3, π/3) ≈ R22(0.3)*sin(2π/3)

        R31(ρ) = 3ρ^3 - 2ρ
        R33(ρ) = ρ^3
        @test TTCal.zernike(3, +1, 0.3, π/3) ≈ R31(0.3)*cos(π/3)
        @test TTCal.zernike(3, +3, 0.3, π/3) ≈ R33(0.3)*cos(π)

        R40(ρ) = 6ρ^4 - 6ρ^2 + 1
        R44(ρ) = ρ^4
        @test TTCal.zernike(4,  0, 0.3, π/3) ≈ R40(0.3)
        @test TTCal.zernike(4, +4, 0.3, π/3) ≈ R44(0.3)*cos(4π/3)

        R60(ρ) = 20ρ^6 - 30ρ^4 + 12ρ^2 - 1
        R64(ρ) = 6ρ^6 - 5ρ^4
        @test TTCal.zernike(6,  0, 0.3, π/3) ≈ R60(0.3)
        @test TTCal.zernike(6, +4, 0.3, π/3) ≈ R64(0.3)*cos(4π/3)
    end
end

