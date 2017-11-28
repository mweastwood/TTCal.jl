@testset "stokes.jl" begin
    @test zero(StokesVector) === StokesVector(0, 0, 0, 0)
    @test one(StokesVector) === StokesVector(1, 0, 0, 0)
    @test repr(StokesVector(1, 2, 3, 4)) == "(1.000, 2.000, 3.000, 4.000)"

    v1 = rand(4)
    v2 = rand(4)
    s1 = StokesVector(v1)
    s2 = StokesVector(v2)
    @test Vector(s1+s2) == v1+v2
    @test Vector(s1-s2) == v1-v2

    a = rand()
    @test Vector(a*s1) == a*v1
    @test Vector(s1*a) == v1*a
    @test Vector(s1/a) == v1/a

    @test TTCal.to_linear*TTCal.to_stokes ≈ eye(4)

    stokes  = rand(StokesVector)
    linear  = HermitianJonesMatrix(stokes)
    stokes′ = StokesVector(linear)
    @test Vector(stokes) ≈ Vector(stokes′)

    stokes = StokesVector(1, 0, 0, 0)
    @test HermitianJonesMatrix(stokes) == HermitianJonesMatrix(1, 0, 1)
    stokes = StokesVector(0, 1, 0, 0)
    @test HermitianJonesMatrix(stokes) == HermitianJonesMatrix(1, 0, -1)
    stokes = StokesVector(0, 0, 1, 0)
    @test HermitianJonesMatrix(stokes) == HermitianJonesMatrix(0, 1, 0)
    stokes = StokesVector(0, 0, 0, 1)
    @test HermitianJonesMatrix(stokes) == HermitianJonesMatrix(0, -1im, 0)

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
    flux   = HermitianJonesMatrix(stokes)
    J = rand(JonesMatrix)
    M = MuellerMatrix(J)

    @test Matrix(TTCal.congruence_transform(J, flux)) ≈ Matrix(J)*Matrix(flux)*Matrix(J')
    @test Matrix(TTCal.congruence_transform(J, flux)) ≈ Matrix(HermitianJonesMatrix(M*stokes))
end

