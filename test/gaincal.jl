let
    N = 100
    data  = rand(Complex64,N,N)
    model = copy(data)

    # Test that the first guess is exact
    # in the ideal situation.
    firstguess = TTCal.gaincal_firstguess(data,model)
    firstguess = firstguess * conj(firstguess[1])/abs(firstguess[1]) # zero the phase
    @test firstguess ≈ ones(Complex64,N)

    # Test that the step is zero when the
    # gains are already at the optimal value.
    gains = ones(Complex64,N)
    @test TTCal.gaincal_step(gains,data,model) ≈ zeros(Complex64,N)
end

let
    Nant = 10
    Nfreq = 20

    cal = GainCalibration(Nant,Nfreq)
    @test TTCal.Nant(cal) == Nant
    @test TTCal.Nfreq(cal) == Nfreq
    @test all(cal.gains .== 1)
    @test all(cal.flags .== false)

    # The inverse of a calibration with unity gains should be itself.
    inverse = TTCal.invert(cal)
    @test cal.gains == inverse.gains
    @test cal.flags == inverse.flags

    # Test that fixphase! actually sets the phase to zero.
    rand!(cal.gains)
    TTCal.fixphase!(cal,1)
    @test maximum(angle(cal.gains[1,:,:])) < eps(Float32)
    TTCal.fixphase!(cal,2)
    @test maximum(angle(cal.gains[2,:,:])) < eps(Float32)
end

function test_gaincal(cal,data,model,
                      ant1,ant2,
                      maxiter,tolerance)
    Nant = TTCal.Nant(cal)
    Nfreq = TTCal.Nfreq(cal)

    # Run as `gaincal!(...)`
    mycal = TTCal.GainCalibration(Nant,Nfreq)
    flags = zeros(Bool,size(data))
    TTCal.gaincal!(mycal,data,model,flags,
                   ant1,ant2,maxiter,tolerance,1)
    @test !any(mycal.flags)
    @test mycal.gains ≈ cal.gains
end

# Unity gains
let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1,ant2 = ant1ant2(Nant)

    cal = TTCal.GainCalibration(Nant,Nfreq)
    data  = Array{Complex64,3}(complex(randn(4,Nfreq,Nbase),randn(4,Nfreq,Nbase)))
    model = copy(data)

    test_gaincal(cal,data,model,ant1,ant2,100,eps(Float32))
end

# Random gains
let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1,ant2 = ant1ant2(Nant)

    cal = TTCal.GainCalibration(Nant,Nfreq)
    for i in eachindex(cal.gains)
        cal.gains[i] = complex(randn(),randn())
    end
    TTCal.fixphase!(cal,1)
    data  = Array{Complex64,3}(complex(randn(4,Nfreq,Nbase),randn(4,Nfreq,Nbase)))
    model = copy(data)
    corrupt!(data,cal,ant1,ant2)

    test_gaincal(cal,data,model,ant1,ant2,100,eps(Float32))
end

# Random gains and corrupted autocorrelations
let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1,ant2 = ant1ant2(Nant)

    cal = TTCal.GainCalibration(Nant,Nfreq)
    for i in eachindex(cal.gains)
        cal.gains[i] = complex(randn(),randn())
    end
    TTCal.fixphase!(cal,1)
    data  = Array{Complex64,3}(complex(randn(4,Nfreq,Nbase),randn(4,Nfreq,Nbase)))
    model = copy(data)
    corrupt!(data,cal,ant1,ant2)
    α = 1
    for ant = 1:Nant
        data[:,:,α] = rand(4,Nfreq)
        α += Nant-ant+1
    end

    test_gaincal(cal,data,model,ant1,ant2,200,eps(Float32))
end

# Test the interface for interacting with measurement sets
let
    Nant = 10
    Nfreq = 2

    # Run as `gaincal(...)`
    name,ms = createms(Nant,Nfreq)
    sources = readsources("sources.json")
    mycal = gaincal(ms,sources,maxiter=100,tolerance=Float64(eps(Float32)))
    @test !any(mycal.flags)
    @test mycal.gains ≈ ones(Complex64,Nant,Nfreq,2)
    unlock(ms)

    # Run from the command line
    output_name = tempname()
    run(`$JULIA_HOME/julia ../src/ttcal.jl gaincal --input $name --output $output_name --sources sources.json --maxiter 100 --tolerance $(eps(Float32))`)
    mycal = TTCal.read(output_name)
    @test !any(mycal.flags)
    @test mycal.gains ≈ ones(Complex64,Nant,Nfreq,2)
end

