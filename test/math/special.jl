@testset "special.jl" begin
    @testset "Hermite polynomials" begin
        H0(x) = 1
        H1(x) = 2x
        H2(x) = 4x^2 - 2
        H3(x) = 8x^3 - 12x
        H4(x) = 16x^4 - 48x^2 + 12
        H5(x) = 32x^5 - 160x^3 + 120x

        for idx = 1:5
            x = randn()
            @test TTCal.hermite(0, x) ≈ H0(x)
            @test TTCal.hermite(1, x) ≈ H1(x)
            @test TTCal.hermite(2, x) ≈ H2(x)
            @test TTCal.hermite(3, x) ≈ H3(x)
            @test TTCal.hermite(4, x) ≈ H4(x)
            @test TTCal.hermite(5, x) ≈ H5(x)
        end
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

