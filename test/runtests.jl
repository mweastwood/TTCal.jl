using TTCal
using Base.Test
using CasaCore.Quanta
using CasaCore.Measures
using CasaCore.Tables

srand(123)
include("jones.jl")
include("sourcemodel.jl")

# Test fringepattern!
ϕ = linspace(1,10,100)
fringe_naive = exp(1im*ϕ)
fringe = TTCal.fringepattern(ϕ[1],ϕ[2]-ϕ[1],length(ϕ))
@test vecnorm(fringe-fringe_naive)/vecnorm(fringe_naive) < 10eps(Float32)

# Define the interferometer
const x = [  0.;  10.;  30.;  70.; 150.]
const y = [150.;  70.;  30.;  10.;   0.]
const z = [  0.;  -1.;  +1.;  -2.;  +2.]
const ν = [40.0e6:10.0e6:60.0e6;]
const t = (2015.-1858.)*365.*24.*60.*60. # a rough, current Julian date (in seconds)

const Nant  = length(x)
const Nbase = div(Nant*(Nant-1),2) + Nant
const Nfreq = length(ν)

const frame = ReferenceFrame()
set!(frame,Epoch(Measures.UTC,Quantity(t,Second)))
set!(frame,observatory("OVRO_MMA"))

const zenith = Direction(Measures.AZEL,Quantity(0.0,Degree),Quantity(90.0,Degree))
const phase_dir = measure(frame,zenith,Measures.J2000)

const pos = observatory("OVRO_MMA")

function xyz2uvw(x,y,z)
    u = Array{Float64}(Nbase)
    v = Array{Float64}(Nbase)
    w = Array{Float64}(Nbase)
    α = 1
    for i = 1:Nant, j = i:Nant
        u[α] = x[j]-x[i]
        v[α] = y[j]-y[i]
        w[α] = z[j]-z[i]
        α += 1
    end
    u,v,w
end
const u,v,w = xyz2uvw(x,y,z)

function getant1ant2()
    ant1 = Array(Int32,Nbase)
    ant2 = Array(Int32,Nbase)
    α = 1
    for i = 1:Nant, j = i:Nant
        ant1[α] = i
        ant2[α] = j
        α += 1
    end
    ant1,ant2
end
const ant1,ant2 = getant1ant2()

function createms()
    name  = tempname()*".ms"
    table = Table(name)

    subtable = Table("$name/SPECTRAL_WINDOW")
    Tables.addRows!(subtable,1)
    subtable["CHAN_FREQ"] = reshape(ν,length(ν),1)
    unlock(subtable)

    subtable = Table("$name/ANTENNA")
    Tables.addRows!(subtable,Nant)
    x,y,z = Measures.xyz_in_meters(pos)
    subtable["POSITION"] = [x;y;z]*ones(1,Nant)
    unlock(subtable)

    subtable = Table("$name/FIELD")
    Tables.addRows!(subtable,1)
    subtable["PHASE_DIR"] = reshape([longitude(phase_dir);latitude(phase_dir)],2,1)
    unlock(subtable)

    Tables.addRows!(table,Nbase)
    table[kw"SPECTRAL_WINDOW"] = "Table: $name/SPECTRAL_WINDOW"
    table[kw"ANTENNA"] = "Table: $name/ANTENNA"
    table[kw"FIELD"] = "Table: $name/FIELD"
    table["ANTENNA1"] = ant1-1
    table["ANTENNA2"] = ant2-1
    table["UVW"] = [u v w]'
    table["TIME"] = fill(float(t),Nbase)

    name,table
end

const maxiter = 100
const tolerance = 1e-4

const gaintable = tempname()*".bcal"
const gaincal_args = Dict("--input"     => "",
                          "--output"    => gaintable,
                          "--sources"   => "sources.json",
                          "--maxiter"   => maxiter,
                          "--tolerance" => tolerance)

const sources = filter(source -> TTCal.isabovehorizon(frame,source),readsources("sources.json"))

function test_gaincal(cal,data,model)
    mycal = TTCal.GainCalibration(Nant,Nfreq)
    flags = zeros(Bool,size(data))
    TTCal.gaincal!(mycal,data,model,flags,
                   ant1,ant2,maxiter,tolerance,1)
    @test vecnorm(mycal.gains-cal.gains)/vecnorm(cal.gains) < 10tolerance

    name,ms = createms()
    gaincal_args["--input"] = name
    ms["DATA"] = data
    ms["FLAG"] = flags
    ms["FLAG_ROW"] = zeros(Bool,Nbase)
    unlock(ms)

    TTCal.run_gaincal(gaincal_args)
    mycal = TTCal.read(gaintable)
    @test vecnorm(mycal.gains-cal.gains)/vecnorm(cal.gains) < 10tolerance

    rm(gaintable)
    run(`$JULIA_HOME/julia ../src/ttcal.jl gaincal --input $name --output $gaintable --sources sources.json --maxiter $maxiter --tolerance $tolerance`)
    mycal = TTCal.read(gaintable)
    @test vecnorm(mycal.gains-cal.gains)/vecnorm(cal.gains) < 10tolerance
end

# Unity gains
function test_one()
    println("1")
    cal = TTCal.GainCalibration(Nant,Nfreq)
    cal.gains[:] = 1
    model = genvis(frame,phase_dir,sources,u,v,w,ν)
    data  = copy(model)
    test_gaincal(cal,data,model)
end
test_one()

# Random gains
function test_two()
    println("2")
    cal = TTCal.GainCalibration(Nant,Nfreq)
    rand!(cal.gains)
    TTCal.fixphase!(cal,1)
    model = genvis(frame,phase_dir,sources,u,v,w,ν)
    data  = copy(model)
    flags = zeros(Bool,size(data))
    corrupt!(data,flags,cal,ant1,ant2)
    test_gaincal(cal,data,model)
end
test_two()

# Random gains
# Corrupted autocorrelations
function test_three()
    println("3")
    cal = TTCal.GainCalibration(Nant,Nfreq)
    rand!(cal.gains)
    TTCal.fixphase!(cal,1)
    model = genvis(frame,phase_dir,sources,u,v,w,ν)
    data  = copy(model)
    flags = zeros(Bool,size(data))
    corrupt!(data,flags,cal,ant1,ant2)
    α = 1
    for ant = 1:Nant
        data[:,:,α] = rand(4,Nfreq)
        α += Nant-ant+1
    end
    test_gaincal(cal,data,model)
end
test_three()

function test_applycal()
    g = 2
    cal = TTCal.GainCalibration(Nant,Nfreq)
    cal.gains[:] = g
    data  = rand(Complex64,4,Nfreq,Nbase)
    flags = zeros(Bool,size(data))

    name,ms = createms()
    ms["DATA"] = data
    ms["FLAG"] = flags
    ms["FLAG_ROW"] = zeros(Bool,Nbase)
    unlock(ms)

    bcal_name = tempname()
    TTCal.write(bcal_name,cal)

    args = Dict("--input" => [name],
                "--calibration" => bcal_name)
    TTCal.run_applycal(args)

    calibrated_ms = Table(name)
    @test calibrated_ms["DATA"] == data / (g*conj(g))
end
test_applycal()

#=
function test_fitvisibilities()
    data = TTCal.visibilities(frame,sources,u,v,w,ν)
    @show sources[1] sources[2]
    for i = 1:length(sources)
        # Perturb each source position by about 0.1 degrees
        dir = sources[i].dir
        ra  = dir.m[1]
        dec = dir.m[2]
        Δra  = 0.1randn()*Degree
        Δdec = 0.1randn()*Degree
        sources[i].dir = Direction(dir.system,(ra+Δra,dec+Δdec))
        # Perturb the flux by about 100 Jy, and the spectral index by about 0.1
        sources[i].flux += 100randn()
        sources[i].index += 0.1randn()
    end
    @show sources[1] sources[2]
    newsources = TTCal.fitvis(frame,data,u,v,w,ν,sources,criteria)
    @show newsources[1] newsources[2]
    newsources = TTCal.fitvis(frame,data,u,v,w,ν,newsources,criteria)
    @show newsources[1] newsources[2]
    newsources = TTCal.fitvis(frame,data,u,v,w,ν,newsources,criteria)
    @show newsources[1] newsources[2]
end
#test_fitvisibilities()
=#

function test_subsrc()
    data  = genvis(frame,phase_dir,sources,u,v,w,ν)
    flags = zeros(Bool,size(data))

    name,ms = createms()
    ms["DATA"] = data
    ms["FLAG"] = flags
    ms["FLAG_ROW"] = zeros(Bool,Nbase)
    subsrc!(ms,sources)

    @test vecnorm(ms["CORRECTED_DATA"]) < eps(Float64)
end
test_subsrc()


