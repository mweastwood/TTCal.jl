using TTCal
using Base.Test
using Base.Dates
using CasaCore.Measures
using CasaCore.Tables
using SIUnits
using NPZ

srand(123)

# Test fringepattern!
ϕ = linspace(1,10,100)
fringe_naive = exp(1im*ϕ)
fringe = TTCal.fringepattern(ϕ[1],ϕ[2]-ϕ[1],length(ϕ))
@test vecnorm(fringe-fringe_naive)/vecnorm(fringe_naive) < 10eps(Float32)

# Define the interferometer
const x = TTCal.addunits([  0.  10.  30.  70. 150.],Meter)
const y = TTCal.addunits([150.  70.  30.  10.   0.],Meter)
const z = TTCal.addunits([  0.  -1.  +1.  -2.  +2.],Meter)
const ν = TTCal.addunits([40.0e6:10.0e6:60.0e6],Hertz)
const t = (2015.-1858.)*365.*24.*60.*60. * Second # a rough, current Julian date

const Nant  = length(x)
const Nbase = div(Nant*(Nant-1),2) + Nant
const Nfreq = length(ν)

function xyz2uvw(x,y,z)
    u = Array(quantity(Float64,Meter),Nbase)
    v = Array(quantity(Float64,Meter),Nbase)
    w = Array(quantity(Float64,Meter),Nbase)
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
    subtable["CHAN_FREQ"] = reshape(TTCal.stripunits(ν),length(ν),1)
    finalize(subtable)

    Tables.addRows!(table,Nbase)
    table[kw"SPECTRAL_WINDOW"] = "Table: $name/SPECTRAL_WINDOW"
    table["ANTENNA1"] = ant1
    table["ANTENNA2"] = ant2
    table["UVW"] = TTCal.stripunits([u v w]')
    table["TIME"] = fill(float(t),Nbase)

    name,table
end

const gaintable = tempname()*".npy"
const args = Dict(      "gaintable" => gaintable,
                  "measurementsets" => [""],
                         "applycal" => false,
                  "doubleprecision" => false,
                            "flags" => Int[],
                            "niter" => 100,
                           "refant" => 1,
                               "RK" => 4,
                              "tol" => 1e-6)

const criteria = TTCal.StoppingCriteria(100,1e-6)
const interferometer = TTCal.Interferometer(Nant,Nfreq,1,Int[])

const frame = ReferenceFrame()
set!(frame,Epoch("UTC",t))
set!(frame,Measures.observatory(frame,"OVRO_MMA"))

const sources = TTCal.getsources(frame)
for source in sources
    @test TTCal.isabovehorizon(frame,source) == true
end

function test_bandpass(gains,data,model)
    mygains = similar(gains)
    TTCal.bandpass!(mygains,data,model,interferometer,criteria)
    @test vecnorm(mygains-gains)/vecnorm(gains) < 1e-5

    name,ms = createms()
    args["measurementsets"] = [name]
    ms["DATA"] = data
    finalize(ms)

    TTCal.run(args)
    mygains = npzread(gaintable)
    @test vecnorm(mygains-gains)/vecnorm(gains) < 1e-5
end

# Unity gains
function test_one()
    println("1")
    gains = ones(Complex64,Nant,2,Nfreq)
    model = TTCal.visibilities(frame,sources,u,v,w,ν)
    data  = copy(model)
    test_bandpass(gains,data,model)
end
test_one()

# Random gains
function test_two()
    println("2")
    gains = rand(Complex64,interferometer.Nant,2,interferometer.Nfreq)
    gains = gains .* conj(gains[1,:,:]) ./ abs(gains[1,:,:])
    model = TTCal.visibilities(frame,sources,u,v,w,ν)
    data  = TTCal.applycal(model,1./gains)
    test_bandpass(gains,data,model)
end
test_two()

# Random gains
# Corrupted autocorrelations
function test_three()
    println("3")
    gains = rand(Complex64,interferometer.Nant,2,interferometer.Nfreq)
    gains = gains .* conj(gains[1,:,:]) ./ abs(gains[1,:,:])
    model = TTCal.visibilities(frame,sources,u,v,w,ν)
    data  = TTCal.applycal(model,1./gains)
    α = 1
    for ant = 1:Nant
        data[:,:,α] = rand(4,Nfreq)
        α += Nant-ant+1
    end
    test_bandpass(gains,data,model)
end
test_three()

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

