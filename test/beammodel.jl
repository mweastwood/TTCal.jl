let
    beam = ConstantBeam()
    for i = 1:3
        ν  = rand() * 100e6
        az = rand() * 2π
        el = rand() * π/2
        @test beam(ν,az,el) == one(JonesMatrix)
    end
end

let α = 2.0
    beam = SineBeam(α)
    for i = 1:3
        ν  = rand() * 100e6
        az = rand() * 2π
        el = rand() * π/2
        M = beam(ν,az,el) |> MuellerMatrix
        @test Matrix(M) ≈ sin(el)^α * eye(4)
    end
end

let
    beam = Memo178Beam()
    @test beam(45e6,0,π/2) == one(JonesMatrix)
end

