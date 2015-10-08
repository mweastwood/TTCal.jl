let
    mat = eye(2,2)
    J = one(JonesMatrix)
    @test Matrix(J) == mat
    @test J == JonesMatrix(mat)
    
    mat = zeros(2,2)
    J = zero(JonesMatrix)
    @test Matrix(J) == mat
    @test J == JonesMatrix(mat)

    mat = rand(Complex64,2,2)
    J = JonesMatrix(mat)
    @test Matrix(J) == mat
    
    mat1 = rand(Complex64,2,2)
    mat2 = rand(Complex64,2,2)
    J1 = JonesMatrix(mat1)
    J2 = JonesMatrix(mat2)
    @test Matrix(J1+J2) == mat1+mat2
    @test Matrix(J1-J2) == mat1-mat2
    @test Matrix(J1*J2) == mat1*mat2
    @test Matrix(J1\J2) ≈ mat1\mat2
    @test Matrix(J1') == mat1'
    @test det(J1) ≈ det(mat1)
    @test Matrix(inv(J1)) ≈ inv(mat1)
    @test norm(J1) ≈ vecnorm(mat1)

    a = rand()
    @test Matrix(a*J1) == a*mat1
    @test Matrix(J1*a) == mat1*a

    @test kron(J1,J2) == kron(mat1,mat2)
end

let
    # test the Mueller matrix generation from the
    # example Jones and Mueller matrices given on:
    # http://scienceworld.wolfram.com/physics/JonesMatrix.html
    # http://scienceworld.wolfram.com/physics/MuellerMatrix.html

    i = 1im

    # linear horizontal polarizer
    J = JonesMatrix(1,0,
                    0,0)
    M = 0.5*[1 1 0 0;
             1 1 0 0;
             0 0 0 0;
             0 0 0 0]
    @test TTCal.mueller(J) ≈ M

    # linear vertical polarizer
    J = JonesMatrix(0,0,
                    0,1)
    M = 0.5*[ 1 -1 0 0;
             -1  1 0 0;
              0  0 0 0;
              0  0 0 0]
    @test TTCal.mueller(J) ≈ M

    # linear polarizer at +45deg
    J = 0.5*JonesMatrix(1,1,
                        1,1)
    M = 0.5*[1 0 1 0;
             0 0 0 0;
             1 0 1 0;
             0 0 0 0]
    @test TTCal.mueller(J) ≈ M

    # linear polarizer at -45deg
    J = 0.5*JonesMatrix( 1,-1,
                        -1, 1)
    M = 0.5*[ 1 0 -1 0;
              0 0  0 0;
             -1 0  1 0;
              0 0  0 0]
    @test TTCal.mueller(J) ≈ M

    # quarter-wave plate, fast axis vertical
    J = exp(i*π/4)*JonesMatrix(1, 0,
                               0,-i)
    M = [1 0 0  0;
         0 1 0  0;
         0 0 0 -1;
         0 0 1  0]
    @test isapprox(TTCal.mueller(J),M,rtol=1e-7)

    # quarter-wave plate, fast axis horizontal
    J = exp(i*π/4)*JonesMatrix(1,0,
                               0,i)
    M = [1 0  0 0;
         0 1  0 0;
         0 0  0 1;
         0 0 -1 0]
    @test isapprox(TTCal.mueller(J),M,rtol=1e-7)

    # circular polarizer, right-handed
    J = 0.5*JonesMatrix( 1,i,
                        -i,1)
    M = 0.5*[1 0 0 1;
             0 0 0 0;
             0 0 0 0;
             1 0 0 1]
    @test TTCal.mueller(J) ≈ M

    # circular polarizer, left-handed
    J = 0.5*JonesMatrix(1,-i,
                        i, 1)
    M = 0.5*[ 1 0 0 -1;
              0 0 0  0;
              0 0 0  0;
             -1 0 0  1]
    @test TTCal.mueller(J) ≈ M
end

