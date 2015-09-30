let
    mat = eye(2,2)
    J = one(JonesMatrix)
    @test convert(Matrix,J) == mat
    @test J == JonesMatrix(mat)
    
    mat = zeros(2,2)
    J = zero(JonesMatrix)
    @test convert(Matrix,J) == mat
    @test J == JonesMatrix(mat)

    mat = rand(Complex64,2,2)
    J = JonesMatrix(mat)
    @test convert(Matrix,J) == mat
    
    mat1 = rand(Complex64,2,2)
    mat2 = rand(Complex64,2,2)
    J1 = JonesMatrix(mat1)
    J2 = JonesMatrix(mat2)
    @test convert(Matrix,J1+J2) == mat1+mat2
    @test convert(Matrix,J1-J2) == mat1-mat2
    @test convert(Matrix,J1*J2) == mat1*mat2
    @test convert(Matrix,J1\J2) ≈ mat1\mat2
    @test convert(Matrix,J1') == mat1'
    @test det(J1) ≈ det(mat1)
    @test convert(Matrix,inv(J1)) ≈ inv(mat1)
    @test norm(J1) ≈ vecnorm(mat1)
end

