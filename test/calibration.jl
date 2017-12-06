@testset "calibration.jl" begin
    Nant  = 256
    Nfreq = 2
    name, ms = createms(Nant, Nfreq)
    meta = TTCal.Metadata(ms)
    beam = TTCal.ConstantBeam()
    sky  = readsky("sources.json")
    #visibilities = genvis(meta, beam, sky)

    #@testset "types" begin
    #    for T in (GainCalibration,PolarizationCalibration)
    #        cal = T(Nant,Nfreq)
    #        el  = eltype(cal.jones)
    #        @test TTCal.Nant(cal) == Nant
    #        @test TTCal.Nfreq(cal) == Nfreq
    #        @test cal.jones == ones(el,Nant,Nfreq)
    #        @test cal.flags == zeros(Bool,Nant,Nfreq)

    #        # Test that the inverse of an identity matrix is itself
    #        rand!(cal.flags)
    #        inverse = TTCal.invert(cal)
    #        @test vecnorm(cal.jones - inverse.jones) == 0
    #        @test cal.flags == inverse.flags
    #    end
    #end

    #@testset "applycal" begin
    #    for T in (GainCalibration,PolarizationCalibration)
    #        @show "applycal", T
    #        g = rand(Complex128)
    #        δ32 = sqrt(eps(Float32))*vecnorm(visibilities.data)
    #        δ64 = sqrt(eps(Float64))*vecnorm(visibilities.data)
    #        random_cal   = T(Nant,Nfreq)
    #        constant_cal = T(Nant,Nfreq)
    #        for i in eachindex(random_cal.jones, constant_cal.jones)
    #            random_cal.jones[i]   =  rand(eltype(  random_cal.jones))
    #            constant_cal.jones[i] = g*one(eltype(constant_cal.jones))
    #        end

    #        # Test that applycal! and corrupt! are inverses
    #        myvis = deepcopy(visibilities)
    #        corrupt!(myvis, meta, random_cal)
    #        applycal!(myvis, meta, random_cal)
    #        @test vecnorm(myvis.data - visibilities.data) < δ64
    #        @test myvis.flags == visibilities.flags

    #        # Test that we get the correct answer with a simple calibration
    #        myvis = deepcopy(visibilities)
    #        applycal!(myvis, meta, constant_cal)
    #        @test vecnorm(myvis.data - visibilities.data/abs2(g)) < δ64
    #        @test myvis.flags == visibilities.flags

    #        # Test that flags are propagated properly
    #        myvis = deepcopy(visibilities)
    #        mycal = deepcopy(constant_cal)
    #        mycal.flags[1,:] = true
    #        expected_flags = zeros(Bool, TTCal.Nbase(meta), Nfreq)
    #        for α = 1:TTCal.Nbase(meta)
    #            antenna1 = meta.baselines[α].antenna1
    #            antenna2 = meta.baselines[α].antenna2
    #            if antenna1 == 1 || antenna2 == 1
    #                expected_flags[α,:] = true
    #            end
    #        end
    #        applycal!(myvis, meta, mycal)
    #        @test vecnorm(myvis.data - visibilities.data/abs2(g)) < δ64
    #        @test myvis.flags == expected_flags

    #        # Test the command line interface
    #        cal_name = tempname()*".jld"
    #        TTCal.write(cal_name, constant_cal)
    #        TTCal.write(ms, "DATA", visibilities)
    #        Tables.unlock(ms)
    #        TTCal.main(["applycal", name, cal_name])
    #        Tables.lock(ms)
    #        myvis = TTCal.read(ms, "DATA")
    #        @test vecnorm(myvis.data - visibilities.data/abs2(g)) < δ32
    #        @test myvis.flags == visibilities.flags
    #    end
    #end

    @testset "stefcal" begin
        for pol in (TTCal.XX, TTCal.YY, TTCal.Dual, TTCal.Full)
            dataset = TTCal.Dataset(meta, polarization=pol)
            rand!(dataset)

            gains = TTCal.Calibration(meta, polarization=pol)
            step  = TTCal.stefcal_step(gains[1, 1].data, dataset[1, 1], dataset[1, 1])
            @test norm(step) < Nant*eps(Float64)
        end
        #N = 100
        #data  = [rand(JonesMatrix) for i = 1:N, j = 1:N]
        #model = copy(data)

        # Test that the step is zero when the Jones matrices are already at the optimal value
        #jones = ones(DiagonalJonesMatrix,N)
        #@test TTCal.stefcal_step(jones,data,model) |> vecnorm < N*eps(Float64)
        #jones = ones(JonesMatrix,N)
        #@test TTCal.stefcal_step(jones,data,model) |> vecnorm < N*eps(Float64)
    end

    #@testset "fixphase" begin
    #    @show "fixphase"
    #    cal = GainCalibration(Nant,Nfreq)
    #    for i in eachindex(cal.jones)
    #        cal.jones[i] = rand(eltype(cal.jones))
    #    end

    #    # Test that fixphase! actually sets the phase to zero.
    #    TTCal.fixphase!(cal,"1x")
    #    for β = 1:Nfreq
    #        @test angle(cal.jones[1,β].xx) < eps(Float64)
    #    end
    #    TTCal.fixphase!(cal,"2y")
    #    for β = 1:Nfreq
    #        @test angle(cal.jones[2,β].yy) < eps(Float64)
    #    end
    #    TTCal.fixphase!(cal,"5y")
    #    for β = 1:Nfreq
    #        @test angle(cal.jones[5,β].yy) < eps(Float64)
    #    end
    #end

    #function test_solve(cal,data,model,maxiter,tolerance)
    #    # Run as `solve!(...)`
    #    println(1)
    #    mycal = similar(cal)
    #    println(1)
    #    TTCal.solve!(mycal,data,model,meta,maxiter,tolerance,true)
    #    println(1)
    #    TTCal.fixphase!(cal,"1x")
    #    println(1)
    #    TTCal.fixphase!(mycal,"1x")
    #    println(1)
    #    jones = cal.jones
    #    println(1)
    #    myjones = mycal.jones
    #    println(1)
    #    jones[cal.flags] = zero(eltype(jones))
    #    println(1)
    #    myjones[mycal.flags] = zero(eltype(myjones))
    #    println(1)
    #    @test cal.flags == mycal.flags
    #    println(1)
    #    @test vecnorm(myjones-jones) < eps(Float32)*vecnorm(jones)
    #end

    #@testset "unity gains" begin
    #    @show "unity gains"
    #    measured = deepcopy(visibilities)
    #    model = deepcopy(visibilities)
    #    for T in (GainCalibration,PolarizationCalibration)
    #        @show T
    #        cal = T(Nant,Nfreq)
    #        test_solve(cal,measured,model,100,eps(Float64))
    #    end
    #end

    #@testset "random gains" begin
    #    @show "random gains"
    #    for T in (GainCalibration,PolarizationCalibration)
    #        @show T
    #        cal = T(Nant,Nfreq)
    #        for i in eachindex(cal.jones)
    #            cal.jones[i] = rand(eltype(cal.jones))
    #        end
    #        measured = Visibilities([rand(JonesMatrix) for α = 1:TTCal.Nbase(meta), β = 1:Nfreq],
    #                                zeros(Bool, TTCal.Nbase(meta), Nfreq))
    #        model = deepcopy(measured)
    #        corrupt!(measured,meta,cal)
    #        test_solve(cal,measured,model,200,10eps(Float64))
    #    end
    #end

    #@testset "corrupted autos" begin
    #    @show "corrupted autos"
    #    for T in (GainCalibration,PolarizationCalibration)
    #        @show T
    #        cal = T(Nant,Nfreq)
    #        for i in eachindex(cal.jones)
    #            cal.jones[i] = rand(eltype(cal.jones))
    #        end
    #        measured = Visibilities([rand(JonesMatrix) for α = 1:TTCal.Nbase(meta), β = 1:Nfreq],
    #                                zeros(Bool, TTCal.Nbase(meta), Nfreq))
    #        model = deepcopy(measured)
    #        corrupt!(measured,meta,cal)
    #        α = 1
    #        for ant = 1:Nant
    #            for β = 1:Nfreq
    #                measured.data[α,β] = rand(JonesMatrix)
    #            end
    #            α += Nant-ant+1
    #        end
    #        test_solve(cal,measured,model,200,10eps(Float64))
    #    end
    #end

    #@testset "flagged antenna" begin
    #    @show "flagged antennas"
    #    for T in (GainCalibration,PolarizationCalibration)
    #        @show T
    #        cal = T(Nant,Nfreq)
    #        for i in eachindex(cal.jones)
    #            cal.jones[i] = rand(eltype(cal.jones))
    #        end
    #        cal.flags[2,:] = true
    #        measured = Visibilities([rand(JonesMatrix) for α = 1:TTCal.Nbase(meta), β = 1:Nfreq],
    #                                zeros(Bool, TTCal.Nbase(meta), Nfreq))
    #        model = deepcopy(measured)
    #        corrupt!(measured,meta,cal)
    #        test_solve(cal,measured,model,200,10eps(Float64))
    #    end
    #end

    #@testset "gaincal" begin
    #    @show "gaincal"
    #    δ = sqrt(eps(Float64))*vecnorm(ones(DiagonalJonesMatrix,Nant,Nfreq))

    #    # Run as `gaincal(...)`
    #    mycal = gaincal(visibilities, meta, beam, sources ,maxiter=100, tolerance=eps(Float64))
    #    TTCal.fixphase!(mycal,"1x")
    #    @test !any(mycal.flags)
    #    @test vecnorm(mycal.jones - ones(DiagonalJonesMatrix,Nant,Nfreq)) < δ

    #    TTCal.write(ms, "DATA", visibilities)
    #    Tables.unlock(ms)

    #    # Run from `main(...)`
    #    output_name = tempname()*".jld"
    #    TTCal.main(["gaincal", name, output_name, "sources.json",
    #                "--maxiter", "100", "--tolerance", "$(eps(Float64))", "--beam","constant"])
    #    mycal = TTCal.read(output_name)
    #    TTCal.fixphase!(mycal,"1x")
    #    @test !any(mycal.flags)
    #    @test vecnorm(mycal.jones - ones(DiagonalJonesMatrix,Nant,Nfreq)) < δ
    #end

    #@testset "polcal" begin
    #    δ = sqrt(eps(Float64))*vecnorm(ones(JonesMatrix,Nant,Nfreq))

    #    # Run as `polcal(...)`
    #    mycal = polcal(visibilities, meta, beam, sources ,maxiter=100, tolerance=eps(Float64))
    #    TTCal.fixphase!(mycal,"1x")
    #    @test !any(mycal.flags)
    #    @test vecnorm(mycal.jones - ones(JonesMatrix,Nant,Nfreq)) < δ

    #    TTCal.write(ms, "DATA", visibilities)
    #    Tables.unlock(ms)

    #    # Run from `main(...)`
    #    output_name = tempname()*".jld"
    #    TTCal.main(["polcal", name, output_name, "sources.json",
    #                "--maxiter", "100", "--tolerance", "$(eps(Float64))", "--beam","constant"])
    #    mycal = TTCal.read(output_name)
    #    TTCal.fixphase!(mycal,"1x")
    #    @test !any(mycal.flags)
    #    @test vecnorm(mycal.jones - ones(JonesMatrix,Nant,Nfreq)) < δ
    #end

#    @testset "numpy I/O" begin
#        let Nfreq = 10, Nant = 5
#            calibration = GainCalibration(Nant, Nfreq)
#            for β = 1:Nfreq, ant = 1:Nant
#                calibration.jones[ant,β] = rand(DiagonalJonesMatrix)
#                calibration.flags[ant,β] = rand(Bool)
#            end
#            filename = tempname()*".npz"
#            TTCal.write_for_python(filename, calibration)
#            test = npzread(filename)
#            test_gains = test["gains"]
#            test_flags = test["flags"]
#            expected_gains = zeros(Complex128, 2, Nant, Nfreq)
#            expected_flags = zeros(Bool, Nant, Nfreq)
#            for β = 1:Nfreq, ant = 1:Nant
#                expected_gains[1, ant, β] = calibration.jones[ant,β].xx
#                expected_gains[2, ant, β] = calibration.jones[ant,β].yy
#                expected_flags[ant, β] = calibration.flags[ant,β]
#            end
#            @test test_gains == expected_gains
#            @test test_flags == expected_flags
#        end
#
#        let Nfreq = 10, Nant = 5
#            calibration = PolarizationCalibration(Nant, Nfreq)
#            for β = 1:Nfreq, ant = 1:Nant
#                calibration.jones[ant,β] = rand(JonesMatrix)
#                calibration.flags[ant,β] = rand(Bool)
#            end
#            filename = tempname()*".npz"
#            TTCal.write_for_python(filename, calibration)
#            test = npzread(filename)
#            test_gains = test["gains"]
#            test_flags = test["flags"]
#            expected_gains = zeros(Complex128, 4, Nant, Nfreq)
#            expected_flags = zeros(Bool, Nant, Nfreq)
#            for β = 1:Nfreq, ant = 1:Nant
#                expected_gains[1, ant, β] = calibration.jones[ant,β].xx
#                expected_gains[2, ant, β] = calibration.jones[ant,β].xy
#                expected_gains[3, ant, β] = calibration.jones[ant,β].yx
#                expected_gains[4, ant, β] = calibration.jones[ant,β].yy
#                expected_flags[ant, β] = calibration.flags[ant,β]
#            end
#            @test test_gains == expected_gains
#            @test test_flags == expected_flags
#        end
#    end

    Tables.delete(ms)
end

