let
    N = 100
    data  = rand(Complex64,N,N)
    model = copy(data)

    # Test that the step is zero when the
    # gain amplitudes are already at the optimal value.
    amps = ones(Float64,N)
    @test TTCal.ampcal_step(amps,data,model) ≈ zeros(Float64,N)
end

let
    Nant = 10
    Nfreq = 20

    cal = AmplitudeCalibration(Nant,Nfreq)
    @test TTCal.Nant(cal) == Nant
    @test TTCal.Nfreq(cal) == Nfreq
    @test all(cal.amplitudes .== 1)
    @test all(cal.flags .== false)

    # The inverse of a calibration with unity amplitudes should be itself.
    rand!(cal.flags)
    inverse = TTCal.invert(cal)
    @test cal.amplitudes == inverse.amplitudes
    @test cal.flags == inverse.flags
end

function test_ampcal(cal,data,model,
                     ant1,ant2,
                     maxiter,tolerance)
    Nant = TTCal.Nant(cal)
    Nfreq = TTCal.Nfreq(cal)

    # Run as `ampcal!(...)`
    mycal = TTCal.AmplitudeCalibration(Nant,Nfreq)
    flags = zeros(Bool,size(data))
    TTCal.ampcal!(mycal,data,model,flags,
                  ant1,ant2,maxiter,tolerance)
    @test !any(mycal.flags)
    @test isapprox(mycal.amplitudes,cal.amplitudes,rtol=sqrt(tolerance))
end

# Unity amplitudes
let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1,ant2 = ant1ant2(Nant)

    cal = TTCal.AmplitudeCalibration(Nant,Nfreq)
    data  = Array{Complex64,3}(complex(randn(4,Nfreq,Nbase),randn(4,Nfreq,Nbase)))
    model = copy(data)

    test_ampcal(cal,data,model,ant1,ant2,100,eps(Float32))
end

# Random amplitudes
let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1,ant2 = ant1ant2(Nant)

    cal = TTCal.AmplitudeCalibration(Nant,Nfreq)
    for i in eachindex(cal.amplitudes)
        cal.amplitudes[i] = rand()
    end
    data  = Array{Complex64,3}(complex(randn(4,Nfreq,Nbase),randn(4,Nfreq,Nbase)))
    model = copy(data)
    corrupt!(data,cal,ant1,ant2)

    test_ampcal(cal,data,model,ant1,ant2,100,eps(Float32))
end

