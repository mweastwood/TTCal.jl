@testset "ionosphere.jl" begin
    # index of refraction is unity when the plasma frequency is zero
    @test TTCal.appleton_hartree(10e6, 0.0) == 1.0
    @test TTCal.appleton_hartree(10e6, 1e6) == sqrt(1 - 0.1^2) < 1
    # the ionosphere is opaque below the plasma frequency
    @test_throws DomainError TTCal.appleton_hartree(1e6, 10e6)

    # test that the ray tracer is finding the point of intersection
    # between lines and circles correctly
    x,y = TTCal.find_intersection_with_circle(0, 0, 1)
    @test x ≈ 1
    @test y ≈ 0
    x,y = TTCal.find_intersection_with_circle(1, 0, 1)
    @test x ≈ 1/sqrt(2)
    @test y ≈ 1/sqrt(2)
    x,y = TTCal.find_intersection_with_circle(1, 0, sqrt(2))
    @test x ≈ 1
    @test y ≈ 1
    x,y = TTCal.find_intersection_with_circle(0, 1/sqrt(2), 1)
    @test x ≈ 1/sqrt(2)
    @test y ≈ 1/sqrt(2)

    # test that the ray tracer against results from Marin's port of Harish's code
    θ  = 0:2:90
    dθ = [-4.21206028e-01,-4.11644598e-01,-3.85108183e-01,-3.46964900e-01,-3.03474035e-01,
          -2.59831790e-01,-2.19380672e-01,-1.83747936e-01,-1.53360973e-01,-1.27960805e-01,
          -1.06972956e-01,-8.97308835e-02,-7.55920944e-02,-6.39883927e-02,-5.44401992e-02,
          -4.65534178e-02,-4.00091364e-02,-3.45514322e-02,-2.99757198e-02,-2.61185950e-02,
          -2.28493966e-02,-2.00633844e-02,-1.76763098e-02,-1.56201328e-02,-1.38396588e-02,
          -1.22899022e-02,-1.09340180e-02,-9.74167760e-03,-8.68778896e-03,-7.75148504e-03,
          -6.91532115e-03,-6.16463541e-03,-5.48703673e-03,-4.87199300e-03,-4.31049805e-03,
          -3.79480094e-03,-3.31818454e-03,-2.87478326e-03,-2.45943189e-03,-2.06753916e-03,
          -1.69498100e-03,-1.33800927e-03,-9.93172742e-04,-6.57247418e-04,-3.27173935e-04,
           0.00000000e+0]
    frequency = 60e6
    plasma_frequency = 6.34887330988e6
    for i = 1:length(θ)
        elevation = deg2rad(θ[i])
        new_elevation = TTCal.ray_trace(elevation, frequency, plasma_frequency)
        answer = deg2rad(dθ[i]) + elevation
        @test new_elevation ≈ answer
    end

    # test the rule of thumb that the plasma frequency is roughly 9√(n)
    # this rule of thumb is good to about 1%
    @test abs(TTCal.plasma_frequency(1) - 9)/9 < 0.01
    @test abs(TTCal.plasma_frequency(100) - 90)/90 < 0.01

    @test TTCal.refract(π/2, frequency, plasma_frequency) == π/2
    for i = 1:length(θ)
        elevation = deg2rad(θ[i])
        new_elevation = TTCal.ray_trace(elevation, frequency, plasma_frequency)
        @test abs(elevation - TTCal.refract(new_elevation, frequency, plasma_frequency)) < deg2rad(0.1/3600)
    end
end

