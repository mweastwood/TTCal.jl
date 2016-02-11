@testset "stokes.jl" begin
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
            @test kron(J1,J2) == kron(mat1,mat2)
        end
    end

    v1 = rand(4)
    v2 = rand(4)
    s1 = StokesVector(v1)
    s2 = StokesVector(v2)
    @test Vector(s1+s2) == v1+v2
    @test Vector(s1-s2) == v1-v2

    a = rand()
    @test Vector(a*s1) == a*v1
    @test Vector(s1*a) == v1*a

    @test TTCal.to_linear*TTCal.to_stokes ≈ eye(4)

    stokes  = rand(StokesVector)
    linear  = TTCal.linear(stokes)
    stokes′ = TTCal.stokes(linear)
    @test Vector(stokes) ≈ Vector(stokes′)

    stokes = StokesVector(1,0,0,0)
    @test TTCal.linear(stokes) == HermitianJonesMatrix(1,0,1)
    stokes = StokesVector(0,1,0,0)
    @test TTCal.linear(stokes) == HermitianJonesMatrix(1,0,-1)
    stokes = StokesVector(0,0,1,0)
    @test TTCal.linear(stokes) == HermitianJonesMatrix(0,1,0)
    stokes = StokesVector(0,0,0,1)
    @test TTCal.linear(stokes) == HermitianJonesMatrix(0,-1im,0)

    # test the Mueller matrix generation from the
    # example Jones and Mueller matrices given on:
    # http://scienceworld.wolfram.com/physics/JonesMatrix.html
    # http://scienceworld.wolfram.com/physics/MuellerMatrix.html

    # note that because I have chosen to apply Jones matrices
    # as JAJ' instead of J'AJ, the sign of some elements of the
    # Jones matrices are swapped

    i = 1im

    # linear horizontal polarizer
    J = JonesMatrix(1,0,
                    0,0)
    M = 0.5*[1 1 0 0;
             1 1 0 0;
             0 0 0 0;
             0 0 0 0] |> MuellerMatrix
    @test norm(MuellerMatrix(J) - M) < eps(Float64)

    # linear vertical polarizer
    J = JonesMatrix(0,0,
                    0,1)
    M = 0.5*[ 1 -1 0 0;
             -1  1 0 0;
              0  0 0 0;
              0  0 0 0] |> MuellerMatrix
    @test norm(MuellerMatrix(J) - M) < eps(Float64)

    # linear polarizer at +45deg
    J = 0.5*JonesMatrix(1,1,
                        1,1)
    M = 0.5*[1 0 1 0;
             0 0 0 0;
             1 0 1 0;
             0 0 0 0] |> MuellerMatrix
    @test norm(MuellerMatrix(J) - M) < eps(Float64)

    # linear polarizer at -45deg
    J = 0.5*JonesMatrix( 1,-1,
                        -1, 1)
    M = 0.5*[ 1 0 -1 0;
              0 0  0 0;
             -1 0  1 0;
              0 0  0 0] |> MuellerMatrix
    @test norm(MuellerMatrix(J) - M) < eps(Float64)

    # quarter-wave plate, fast axis vertical
    J = exp(i*π/4)*JonesMatrix(1,0,
                               0,i)
    M = [1 0 0  0;
         0 1 0  0;
         0 0 0 -1;
         0 0 1  0] |> MuellerMatrix
    @test norm(MuellerMatrix(J) - M) < eps(Float64)

    # quarter-wave plate, fast axis horizontal
    J = exp(i*π/4)*JonesMatrix(1, 0,
                               0,-i)
    M = [1 0  0 0;
         0 1  0 0;
         0 0  0 1;
         0 0 -1 0] |> MuellerMatrix
    @test norm(MuellerMatrix(J) - M) < eps(Float64)

    # circular polarizer, right-handed
    J = 0.5*JonesMatrix(1,-i,
                        i, 1)
    M = 0.5*[1 0 0 1;
             0 0 0 0;
             0 0 0 0;
             1 0 0 1] |> MuellerMatrix
    @test norm(MuellerMatrix(J) - M) < eps(Float64)

    # circular polarizer, left-handed
    J = 0.5*JonesMatrix( 1,i,
                        -i,1)
    M = 0.5*[ 1 0 0 -1;
              0 0 0  0;
              0 0 0  0;
             -1 0 0  1] |> MuellerMatrix
    @test norm(MuellerMatrix(J) - M) < eps(Float64)

    stokes = rand(StokesVector)
    flux   = TTCal.linear(stokes)
    J = rand(JonesMatrix)
    M = MuellerMatrix(J)

    @test Matrix(TTCal.congruence_transform(J,flux)) ≈ Matrix(J)*Matrix(flux)*Matrix(J')
    @test Matrix(TTCal.congruence_transform(J,flux)) ≈ Matrix(TTCal.linear(M*stokes))
end

