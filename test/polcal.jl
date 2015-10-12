let
    N = 100
    data  = [rand(JonesMatrix) for i = 1:N, j = 1:N]
    model = copy(data)

    # Test that the step is zero when the
    # Jones matrices are already at the optimal value.
    jones = ones(JonesMatrix,N)
    @test norm(TTCal.polcal_step(jones,data,model),Inf) < 5eps(Float32)
end

let
    Nant = 10
    Nfreq = 20

    cal = PolarizationCalibration(Nant,Nfreq)
    @test TTCal.Nant(cal) == Nant
    @test TTCal.Nfreq(cal) == Nfreq
    @test all(cal.jones .== one(JonesMatrix))
    @test all(cal.flags .== false)

    # The inverse of a calibration with identity Jones matrices should be itself.
    inverse = TTCal.invert(cal)
    @test vecnorm(cal.jones-inverse.jones,Inf) < 5eps(Float32)
    @test cal.flags == inverse.flags
end

