let
    N = 100
    data  = rand(Complex64,N,N)
    model = copy(data)

    # Test that the step is zero when the
    # gains are already at the optimal value.
    gains = ones(Complex64,N)
    @test TTCal.gaincal_step(gains,data,model) ≈ zeros(Complex64,N)

    amps = ones(Float64,N)
    @test TTCal.ampcal_step(amps,data,model) ≈ zeros(Float64,N)
end

let
    Nant = 10
    Nfreq = 20

    for T in (GainCalibration,AmplitudeCalibration)
        cal = T(Nant,Nfreq)
        @test TTCal.Nant(cal) == Nant
        @test TTCal.Nfreq(cal) == Nfreq
        @test all(cal.gains .== 1)
        @test all(cal.flags .== false)

        # The inverse of a calibration with unity gains should be itself.
        rand!(cal.flags)
        inverse = TTCal.invert(cal)
        @test cal.gains == inverse.gains
        @test cal.flags == inverse.flags
    end
end

let
    Nant = 10
    Nfreq = 20
    cal = GainCalibration(Nant,Nfreq)

    # Test that fixphase! actually sets the phase to zero.
    rand!(cal.gains)
    TTCal.fixphase!(cal,"1x")
    @test maximum(angle(cal.gains[1,1,:])) < eps(Float32)
    TTCal.fixphase!(cal,"2y")
    @test maximum(angle(cal.gains[2,2,:])) < eps(Float32)
    TTCal.fixphase!(cal,"5y")
    @test maximum(angle(cal.gains[2,5,:])) < eps(Float32)
end

function test_solve(cal,data,model,
                    ant1,ant2,maxiter,tolerance)
    Nant = TTCal.Nant(cal)
    Nfreq = TTCal.Nfreq(cal)

    # Run as `solve!(...)`
    mycal = similar(cal)
    flags = zeros(Bool,size(data))
    TTCal.solve!(mycal,data,model,flags,
                 ant1,ant2,maxiter,tolerance,"1x")
    @test !any(mycal.flags)
    @test isapprox(mycal.gains,cal.gains,atol=sqrt(tolerance))
end

# Unity gains
let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1,ant2 = ant1ant2(Nant)

    data  = Array{Complex64,3}(complex(randn(4,Nfreq,Nbase),randn(4,Nfreq,Nbase)))
    model = copy(data)

    for T in (GainCalibration,AmplitudeCalibration)
        cal = T(Nant,Nfreq)
        test_solve(cal,data,model,ant1,ant2,100,eps(Float32))
    end
end

# Random gains
let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1,ant2 = ant1ant2(Nant)

    for T in (GainCalibration,AmplitudeCalibration)
        cal = T(Nant,Nfreq)
        rand!(cal.gains)
        TTCal.fixphase!(cal,"1x")
        data  = Array{Complex64,3}(complex(randn(4,Nfreq,Nbase),randn(4,Nfreq,Nbase)))
        model = copy(data)
        corrupt!(data,cal,ant1,ant2)
        test_solve(cal,data,model,ant1,ant2,100,eps(Float32))
    end
end

# Random gains and corrupted autocorrelations
let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1,ant2 = ant1ant2(Nant)

    for T in (GainCalibration,AmplitudeCalibration)
        cal = T(Nant,Nfreq)
        rand!(cal.gains)
        TTCal.fixphase!(cal,"1x")
        data  = Array{Complex64,3}(complex(randn(4,Nfreq,Nbase),randn(4,Nfreq,Nbase)))
        model = copy(data)
        corrupt!(data,cal,ant1,ant2)
        α = 1
        for ant = 1:Nant
            data[:,:,α] = rand(4,Nfreq)
            α += Nant-ant+1
        end
        test_solve(cal,data,model,ant1,ant2,200,eps(Float32))
    end

end

# Test the interface for interacting with measurement sets
let
    Nant = 10
    Nfreq = 2

    # Run as `gaincal(...)`
    name,ms = createms(Nant,Nfreq)
    sources = readsources("sources.json")
    mycal = gaincal(ms,sources,TTCal.ConstantBeam(),
                    maxiter=100,tolerance=Float64(eps(Float32)))
    @test !any(mycal.flags)
    @test mycal.gains ≈ ones(Complex64,2,Nant,Nfreq)
    unlock(ms)

    # Run from `main(...)`
    output_name = tempname()
    TTCal.main(["gaincal","--input",name,"--output",output_name,"--sources","sources.json","--maxiter","100","--tolerance","$(eps(Float32))"])
    mycal = TTCal.read(output_name)
    @test !any(mycal.flags)
    @test mycal.gains ≈ ones(Complex64,2,Nant,Nfreq)
end

