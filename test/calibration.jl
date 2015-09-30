# Test that applycal! and corrupt! are inverses
let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1,ant2 = ant1ant2(Nant)

    data  = rand(Complex64,4,Nfreq,Nbase)
    data′ = copy(data)
    flags = zeros(Bool,4,Nfreq,Nbase)

    # ampcal
    cal = TTCal.AmplitudeCalibration(Nant,Nfreq)
    for i in eachindex(cal.amplitudes)
        cal.amplitudes[i] = rand()
    end
    corrupt!(data,flags,cal,ant1,ant2)
    applycal!(data,flags,cal,ant1,ant2)
    @test data ≈ data′

    # gaincal
    cal = TTCal.GainCalibration(Nant,Nfreq)
    for i in eachindex(cal.gains)
        cal.gains[i] = complex(randn(),randn())
    end
    corrupt!(data,flags,cal,ant1,ant2)
    applycal!(data,flags,cal,ant1,ant2)
    @test data ≈ data′
end

# Test the interface for interacting with measurement sets
let
    Nant = 10
    Nfreq = 2

    g = 2
    cal = TTCal.GainCalibration(Nant,Nfreq)
    cal.gains[:] = g

    # Run as `applycal(...)`
    name,ms = createms(Nant,Nfreq)
    data  = ms["DATA"]
    applycal!(ms,cal)
    data′ = ms["DATA"]
    @test data/(g*conj(g)) ≈ data′
    unlock(ms)

    # Run from the command line
    cal_name = tempname()
    TTCal.write(cal_name,cal)
    run(`$JULIA_HOME/julia ../src/ttcal.jl applycal --input $name --calibration $cal_name`)
    calibrated_ms = Table(name)
    @test data/(g*conj(g)) ≈ calibrated_ms["DATA"]
end

