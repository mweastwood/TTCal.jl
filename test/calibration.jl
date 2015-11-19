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

let
    N = 100
    data  = [rand(JonesMatrix) for i = 1:N, j = 1:N]
    model = copy(data)

    # Test that the step is zero when the
    # Jones matrices are already at the optimal value.
    jones = ones(DiagonalJonesMatrix,N)
    @test TTCal.stefcal_step(jones,data,model) |> vecnorm < N*eps(Float64)
    jones = ones(JonesMatrix,N)
    @test TTCal.stefcal_step(jones,data,model) |> vecnorm < N*eps(Float64)
end

let
    Nant = 10
    Nfreq = 20

    for T in (GainCalibration,PolarizationCalibration)
        cal = T(Nant,Nfreq)
        el  = eltype(cal.jones)
        @test TTCal.Nant(cal) == Nant
        @test TTCal.Nfreq(cal) == Nfreq
        @test cal.jones == ones(el,Nant,Nfreq)
        @test cal.flags == zeros(Bool,Nant,Nfreq)

        # The inverse of a calibration with identity Jones matrices should be itself.
        rand!(cal.flags)
        inverse = TTCal.invert(cal)
        @test vecnorm(cal.jones - inverse.jones) == 0
        @test cal.flags == inverse.flags
    end
end

let
    Nant = 10
    Nfreq = 3
    cal = GainCalibration(Nant,Nfreq)
    for i in eachindex(cal.jones)
        cal.jones[i] = rand(eltype(cal.jones))
    end

    # Test that fixphase! actually sets the phase to zero.

    TTCal.fixphase!(cal,"1x")
    for β = 1:Nfreq
        @test angle(cal.jones[1,β].xx) < eps(Float64)
    end

    TTCal.fixphase!(cal,"2y")
    for β = 1:Nfreq
        @test angle(cal.jones[2,β].yy) < eps(Float64)
    end

    TTCal.fixphase!(cal,"5y")
    for β = 1:Nfreq
        @test angle(cal.jones[5,β].yy) < eps(Float64)
    end
end

function test_solve(cal,data,model,
                    ant1,ant2,maxiter,tolerance)
    Nant = TTCal.Nant(cal)
    Nfreq = TTCal.Nfreq(cal)

    # Run as `solve!(...)`
    mycal = similar(cal)
    flags = zeros(Bool,size(data))
    TTCal.solve!(mycal,data,model,flags,
                 ant1,ant2,maxiter,tolerance)
    TTCal.fixphase!(mycal,"1x")
    @test !any(mycal.flags)
    @test vecnorm(mycal.jones-cal.jones) < eps(Float32)
end

# Unity gains
let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1,ant2 = ant1ant2(Nant)

    data  = Array{Complex64,3}(complex(randn(4,Nfreq,Nbase),randn(4,Nfreq,Nbase)))
    model = copy(data)

    for T in (GainCalibration,PolarizationCalibration)
        cal = T(Nant,Nfreq)
        test_solve(cal,data,model,ant1,ant2,100,eps(Float64))
    end
end

# Random gains
let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1,ant2 = ant1ant2(Nant)

    for T in (GainCalibration,PolarizationCalibration)
        cal = T(Nant,Nfreq)
        for i in eachindex(cal.jones)
            cal.jones[i] = rand(eltype(cal.jones))
        end
        TTCal.fixphase!(cal,"1x")
        data  = Array{Complex64,3}(complex(randn(4,Nfreq,Nbase),randn(4,Nfreq,Nbase)))
        model = copy(data)
        corrupt!(data,cal,ant1,ant2)
        test_solve(cal,data,model,ant1,ant2,200,eps(Float64))
    end
end

# Random gains and corrupted autocorrelations
let
    Nant = 10
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1,ant2 = ant1ant2(Nant)

    for T in (GainCalibration,PolarizationCalibration)
        cal = T(Nant,Nfreq)
        for i in eachindex(cal.jones)
            cal.jones[i] = rand(eltype(cal.jones))
        end
        TTCal.fixphase!(cal,"1x")
        data  = Array{Complex64,3}(complex(randn(4,Nfreq,Nbase),randn(4,Nfreq,Nbase)))
        model = copy(data)
        corrupt!(data,cal,ant1,ant2)
        α = 1
        for ant = 1:Nant
            data[:,:,α] = rand(4,Nfreq)
            α += Nant-ant+1
        end
        test_solve(cal,data,model,ant1,ant2,200,eps(Float64))
    end

end

# Test the interface for interacting with measurement sets
let
    Nant = 10
    Nfreq = 2

    # Run as `gaincal(...)`
    name,ms = createms(Nant,Nfreq)
    sources = readsources("sources.json")
    ms.table["DATA"] = genvis(ms,sources,ConstantBeam())

    mycal = gaincal(ms,sources,ConstantBeam(),
                    maxiter=100,tolerance=eps(Float64))
    @test !any(mycal.flags)
    @test vecnorm(mycal.jones - ones(DiagonalJonesMatrix,Nant,Nfreq)) < sqrt(eps(Float64))
    unlock(ms)

    # Run from `main(...)`
    output_name = tempname()*".jld"
    TTCal.main(["gaincal","--input",name,"--output",output_name,"--sources","sources.json","--maxiter","100","--tolerance","$(eps(Float64))","--beam","constant"])
    mycal = TTCal.read(output_name)
    @test !any(mycal.flags)
    @test vecnorm(mycal.jones - ones(DiagonalJonesMatrix,Nant,Nfreq)) < sqrt(eps(Float64))
end

