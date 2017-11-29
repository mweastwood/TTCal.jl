@testset "jones.jl" begin
    for T in (TTCal.JonesMatrix, TTCal.DiagonalJonesMatrix, TTCal.HermitianJonesMatrix)
        @test zero(T) |> Matrix == [0 0; 0 0]
        @test  one(T) |> Matrix == [1 0; 0 1]

        a = rand(Complex128)
        J1 = rand(T)
        J2 = rand(T)
        mat1 = Matrix(J1)
        mat2 = Matrix(J2)
        @test a*J1 == a*mat1
        @test J1*a == mat1*a
        @test J1/a ≈ mat1/a
        @test J1./a ≈ mat1/a
        @test J1+J2 == mat1+mat2
        @test J1-J2 == mat1-mat2
        @test J1*J2 == mat1*mat2
        @test J1\J2 ≈ mat1\mat2
        @test J1/J2 ≈ mat1/mat2
        @test det(J1) ≈ det(mat1)
        @test inv(J1) ≈ inv(mat1)
        @test norm(J1) ≈ norm(mat1)
        @test J1' == mat1'
        @test conj(J1) == conj(mat1)

        #if T == TTCal.JonesMatrix
        #    @test kron(J1, J2) == kron(mat1, mat2)
        #end
    end

    J1 = rand(TTCal.HermitianJonesMatrix)
    J2 = rand(TTCal.JonesMatrix)
    mat1 = Matrix(J1)
    mat2 = Matrix(J2)
    @test J1*J2 == mat1*mat2
    @test J2*J1 == mat2*mat1
    #@test TTCal.congruence_transform(J2, J1) ≈ mat2*mat1*mat2'
    #@test TTCal.make_hermitian(J2) == 0.5*(mat2+mat2')
end

