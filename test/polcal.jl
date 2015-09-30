let
    N = 100
    data  = [JonesMatrix(rand(Complex64,2,2)) for i = 1:N, j = 1:N]
    model = copy(data)

    # Test that the step is zero when the
    # Jones matrices are already at the optimal value.
    jones = ones(JonesMatrix,N)
    @test norm(TTCal.polcal_step(jones,data,model),Inf) < 5eps(Float32)
end

