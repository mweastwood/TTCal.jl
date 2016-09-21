@testset "jones.jl" begin
    for T in (JonesMatrix, DiagonalJonesMatrix, HermitianJonesMatrix)
        @test zero(T) |> Matrix == [0 0; 0 0]
        @test  one(T) |> Matrix == [1 0; 0 1]

        a = rand(Complex128)
        J1 = rand(T)
        J2 = rand(T)
        mat1 = Matrix(J1)
        mat2 = Matrix(J2)
        @test Matrix(a*J1) == a*mat1
        @test Matrix(J1*a) == mat1*a
        @test Matrix(J1/a) ≈ mat1/a
        @test Matrix(J1./a) ≈ mat1/a
        @test Matrix(J1+J2) == mat1+mat2
        @test Matrix(J1-J2) == mat1-mat2
        @test Matrix(J1*J2) == mat1*mat2
        @test Matrix(J1\J2) ≈ mat1\mat2
        @test Matrix(J1/J2) ≈ mat1/mat2
        @test det(J1) ≈ det(mat1)
        @test Matrix(inv(J1)) ≈ inv(mat1)
        @test norm(J1) ≈ vecnorm(mat1)
        @test Matrix(J1') == mat1'
        @test Matrix(conj(J1)) == conj(mat1)

        if T == JonesMatrix
            @test kron(J1, J2) == kron(mat1, mat2)
        end
    end

    J1 = rand(HermitianJonesMatrix)
    J2 = rand(JonesMatrix)
    mat1 = Matrix(J1)
    mat2 = Matrix(J2)
    @test Matrix(J1*J2) == mat1*mat2
    @test Matrix(J2*J1) == mat2*mat1
    @test Matrix(TTCal.congruence_transform(J2, J1)) ≈ mat2*mat1*mat2'
    @test Matrix(TTCal.make_hermitian(J2)) == 0.5*(mat2+mat2')
end

