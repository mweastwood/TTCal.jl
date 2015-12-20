@testset "Calibration Tests" begin
    Nant  = 100
    Nfreq = 2
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1,ant2 = ant1ant2(Nant)

    @testset "types" begin
        for T in (GainCalibration,PolarizationCalibration)
            cal = T(Nant,Nfreq)
            el  = eltype(cal.jones)
            @test TTCal.Nant(cal) == Nant
            @test TTCal.Nfreq(cal) == Nfreq
            @test cal.jones == ones(el,Nant,Nfreq)
            @test cal.flags == zeros(Bool,Nant,Nfreq)

            # Test that the inverse of an identity matrix is itself
            rand!(cal.flags)
            inverse = TTCal.invert(cal)
            @test vecnorm(cal.jones - inverse.jones) == 0
            @test cal.flags == inverse.flags
        end
    end

    @testset "applycal" begin
        data  = rand(Complex64,4,Nfreq,Nbase)
        data′ = copy(data)
        flags = zeros(Bool,4,Nfreq,Nbase)

        for T in (GainCalibration,PolarizationCalibration)
            g = rand(Complex64)
            random_cal   = T(Nant,Nfreq)
            constant_cal = T(Nant,Nfreq)
            for i in eachindex(random_cal.jones, constant_cal.jones)
                random_cal.jones[i]   =  rand(eltype(  random_cal.jones))
                constant_cal.jones[i] = g*one(eltype(constant_cal.jones))
            end

            # Test that applycal! and corrupt! are inverses
            data = copy(data′)
            corrupt!(data,flags,random_cal,ant1,ant2)
            applycal!(data,flags,random_cal,ant1,ant2)
            @test data ≈ data′

            # Test that we get the correct answer with a simple calibration
            data = copy(data′)
            applycal!(data,flags,constant_cal,ant1,ant2)
            @test data ≈ data′/(g*conj(g))

            # Test the interface for interacting with measurement sets
            data = copy(data′)
            name,ms = createms(Nant,Nfreq)
            ms.table["DATA"] = data
            applycal!(ms,constant_cal)
            data = TTCal.get_data(ms)
            @test data ≈ data′/(g*conj(g))

            # Test the command line interface
            data = copy(data′)
            cal_name = tempname()*".jld"
            TTCal.write(cal_name,constant_cal)
            name,ms = createms(Nant,Nfreq)
            ms.table["DATA"] = data
            TTCal.main(["applycal","--input",name,"--calibration",cal_name])
            data = TTCal.get_data(ms)
            @test data ≈ data′/(g*conj(g))
        end
    end

    @testset "stefcal" begin
        N = 100
        data  = [rand(JonesMatrix) for i = 1:N, j = 1:N]
        model = copy(data)

        # Test that the step is zero when the Jones matrices are already at the optimal value
        jones = ones(DiagonalJonesMatrix,N)
        @test TTCal.stefcal_step(jones,data,model) |> vecnorm < N*eps(Float64)
        jones = ones(JonesMatrix,N)
        @test TTCal.stefcal_step(jones,data,model) |> vecnorm < N*eps(Float64)
    end

    @testset "fixphase" begin
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

    function test_solve(cal,data,model,maxiter,tolerance)
        # Run as `solve!(...)`
        mycal = similar(cal)
        flags = zeros(Bool,size(data))
        TTCal.solve!(mycal,data,model,flags,ant1,ant2,maxiter,tolerance,false,quiet=true)
        TTCal.fixphase!(mycal,"1x")
        @test !any(mycal.flags)
        @test vecnorm(mycal.jones-cal.jones) < eps(Float32)*vecnorm(cal.jones)
    end

    @testset "unity gains" begin
        data  = Array{Complex64,3}(complex(randn(4,Nfreq,Nbase),randn(4,Nfreq,Nbase)))
        model = copy(data)
        for T in (GainCalibration,PolarizationCalibration)
            cal = T(Nant,Nfreq)
            test_solve(cal,data,model,100,eps(Float64))
        end
    end

    @testset "random gains" begin
        for T in (GainCalibration,PolarizationCalibration)
            cal = T(Nant,Nfreq)
            for i in eachindex(cal.jones)
                cal.jones[i] = rand(eltype(cal.jones))
            end
            TTCal.fixphase!(cal,"1x")
            data  = Array{Complex64,3}(complex(randn(4,Nfreq,Nbase),randn(4,Nfreq,Nbase)))
            model = copy(data)
            corrupt!(data,cal,ant1,ant2)
            test_solve(cal,data,model,200,eps(Float64))
        end
    end

    @testset "corrupted autos" begin
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
            test_solve(cal,data,model,200,eps(Float64))
        end
    end

    @testset "one bad antenna" begin
        data  = Array{Complex64,3}(complex(randn(4,Nfreq,Nbase),randn(4,Nfreq,Nbase)))
        model = copy(data)
        for α = 1:Nant
            data[:,:,α] = 0
        end
        for T in (GainCalibration,PolarizationCalibration)
            cal = T(Nant,Nfreq)
            mycal = T(Nant,Nfreq)
            flags = zeros(Bool,size(data))
            TTCal.solve!(mycal,data,model,flags,ant1,ant2,200,eps(Float64),true)
            TTCal.fixphase!(mycal,"1y")
            @test all(mycal.flags[1,:])
            @test !any(mycal.flags[2:end,:])
            @test vecnorm(mycal.jones[2:end,:]-cal.jones[2:end,:]) < eps(Float32)*vecnorm(cal.jones[2:end,:])
        end
    end

    @testset "gaincal" begin
        # Run as `gaincal(...)`
        name,ms = createms(Nant,Nfreq)
        sources = readsources("sources.json")
        ms.table["DATA"] = genvis(ms,sources,ConstantBeam())

        mycal = gaincal(ms,sources,ConstantBeam(),maxiter=100,tolerance=eps(Float64))
        @test !any(mycal.flags)
        @test vecnorm(mycal.jones - ones(DiagonalJonesMatrix,Nant,Nfreq)) < sqrt(eps(Float64))
        unlock(ms)

        # Run from `main(...)`
        output_name = tempname()*".jld"
        TTCal.main(["gaincal","--input",name,"--output",output_name,"--sources","sources.json",
                    "--maxiter","100","--tolerance","$(eps(Float64))","--beam","constant"])
        mycal = TTCal.read(output_name)
        @test !any(mycal.flags)
        @test vecnorm(mycal.jones - ones(DiagonalJonesMatrix,Nant,Nfreq)) < sqrt(eps(Float64))
    end

    @testset "polcal" begin
        # Run as `polcal(...)`
        name,ms = createms(Nant,Nfreq)
        sources = readsources("sources.json")
        ms.table["DATA"] = genvis(ms,sources,ConstantBeam())

        mycal = polcal(ms,sources,ConstantBeam(),maxiter=100,tolerance=eps(Float64))
        @test !any(mycal.flags)
        @test vecnorm(mycal.jones - ones(JonesMatrix,Nant,Nfreq)) < sqrt(eps(Float64))
        unlock(ms)

        # Run from `main(...)`
        output_name = tempname()*".jld"
        TTCal.main(["polcal","--input",name,"--output",output_name,"--sources","sources.json",
                    "--maxiter","100","--tolerance","$(eps(Float64))","--beam","constant"])
        mycal = TTCal.read(output_name)
        @test !any(mycal.flags)
        @test vecnorm(mycal.jones - ones(JonesMatrix,Nant,Nfreq)) < sqrt(eps(Float64))
    end
end

