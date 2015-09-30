let
    N = 100
    data  = rand(Complex64,N,N)
    model = copy(data)

    # Test that the step is zero when the
    # gain amplitudes are already at the optimal value.
    amps = ones(Float64,N)
    @test TTCal.ampcal_step(amps,data,model) ≈ zeros(Float64,N)
end

