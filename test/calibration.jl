# Test that applycal! and corrupt! are inverses
let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1,ant2 = ant1ant2(Nant)

    data  = rand(Complex64,4,Nfreq,Nbase)
    data′ = copy(data)
    flags = zeros(Bool,4,Nfreq,Nbase)

    for T in (GainCalibration,PolarizationCalibration)
        cal = T(Nant,Nfreq)
        for i in eachindex(cal.jones)
            cal.jones[i] = rand(eltype(cal.jones))
        end
        corrupt!(data,flags,cal,ant1,ant2)
        applycal!(data,flags,cal,ant1,ant2)
        @test data ≈ data′
    end
end

# Test the interface for interacting with measurement sets
let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant

    g = 2
    cal = GainCalibration(Nant,Nfreq)
    for i in eachindex(cal.jones)
        cal.jones[i] = DiagonalJonesMatrix(g,g)
    end

    # Run as `applycal(...)`
    name,ms = createms(Nant,Nfreq)
    ms.table["DATA"] = rand(Complex64,4,Nfreq,Nbase)
    data  = TTCal.get_data(ms)
    applycal!(ms,cal)
    data′ = TTCal.get_data(ms)
    @test data/(g*conj(g)) ≈ data′
    unlock(ms)

    # Run from `main(...)`
    cal_name = tempname()*".jld"
    TTCal.write(cal_name,cal)
    TTCal.main(["applycal","--input",name,"--calibration",cal_name])
    calibrated_ms = Table(name)
    @test data′/(g*conj(g)) ≈ calibrated_ms["DATA"]
end

