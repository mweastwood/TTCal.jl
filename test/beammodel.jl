@testset "beammodel.jl" begin
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

    @testset "beam models" begin
        beam = ConstantBeam()
        for i = 1:3
            ν  = rand() * 100e6
            az = rand() * 2π
            el = rand() * π/2
            @test beam(ν,az,el) == one(JonesMatrix)
        end

        α = 2.0
        beam = SineBeam(α)
        for i = 1:3
            ν  = rand() * 100e6
            az = rand() * 2π
            el = rand() * π/2
            M = beam(ν,az,el) |> MuellerMatrix
            @test Matrix(M) ≈ sin(el)^α * eye(4)
        end

        beam = Memo178Beam()
        @test beam(45e6,0,π/2) == one(JonesMatrix)
    end
end

