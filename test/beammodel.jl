@testset "Beam Tests" begin
    @testset "dipole projection" begin
        @test TTCal.dipole_projection(0,π/2) |> Matrix ≈ [1 0; 0 1]
        @test TTCal.dipole_projection(0,0) |> Matrix ≈ [0 0; 0 1]
        @test TTCal.dipole_projection(π/2,0) |> Matrix ≈ [1 0; 0 0]
        @test TTCal.dipole_projection(π,0) |> Matrix ≈ [0 0; 0 1]
        @test TTCal.dipole_projection(3π/2,0) |> Matrix ≈ [1 0; 0 0]
    end

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
    end

    @testset "memo 178 beam" begin
        beam = Memo178Beam()
        @test beam(45e6,0,π/2) == one(JonesMatrix)
    end
end

