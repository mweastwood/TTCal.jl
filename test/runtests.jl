using TTCal
using Base.Test
using Base.Dates
using CasaCore.Tables
using NPZ

function create_ms(name,x,y,z,ν)
    Nant  = length(x)
    Nfreq = length(ν)
    Nbase = div(Nant*(Nant-1),2) + Nant

    ant1 = Array(Cint,Nbase)
    ant2 = Array(Cint,Nbase)
    u = Array(Cdouble,Nbase)
    v = Array(Cdouble,Nbase)
    w = Array(Cdouble,Nbase)

    α = 1
    for i = 1:Nant, j = i:Nant
        ant1[α]  = i
        ant2[α]  = j
        u[α] = x[j]-x[i]
        v[α] = y[j]-y[i]
        w[α] = z[j]-z[i]
        α += 1
    end

    t = (2015.-1858.)*365.*24.*60.*60. # a rough, current Julian date

    table = Table(name)
    subtable = Table("$name/SPECTRAL_WINDOW")

    Tables.addRows!(subtable,1)
    subtable["CHAN_FREQ"] = reshape(ν,length(ν),1)
    finalize(subtable)

    Tables.addRows!(table,Nbase)
    table[kw"SPECTRAL_WINDOW"] = "Table: $name/SPECTRAL_WINDOW"
    table["ANTENNA1"] = ant1
    table["ANTENNA2"] = ant2
    table["UVW"] = [u v w]'
    table["TIME"] = fill(t,Nbase)

    table
end

function test_one()
    # Define the interferometer
    # (we'll assume the antennas all lie on east-west baselines)
    x = [  0.  10.  30.  70. 150.] # m
    y = [  0.   0.   0.   0.   0.] # m
    z = [  0.   0.   0.   0.   0.] # m
    ν = linspace(40e6,60e6,3) # Hz

    name = tempname()*".ms"
    @show name
    ms = create_ms(name,x,y,z,ν)
    interferometer = TTCal.Interferometer(length(x),length(ν),1,Int[])
    sources = TTCal.getsources(ms)
    data = TTCal.visibilities(interferometer,ms,sources)
    ms["DATA"] = data

    gains = TTCal.bandpass(interferometer,[ms],sources,TTCal.BandpassOptions(30,1e-5,4,true))
    truegains = ones(size(gains))

    # TODO: tighten this tolerance (at the moment I'm setting
    # it based on the tolerance passed to BandpassOptions)
    @test vecnorm(gains-truegains)/vecnorm(truegains) < 1e-5

    finalize(ms)

    gaintable = tempname()*".npy"
    args = Dict("gaintable"=>gaintable,
                "measurementsets"=>[name],
                "applycal"=>false,
                "doubleprecision"=>false,
                "flags"=>Int[],
                "niter"=>30,
                "refant"=>1,
                "RK"=>4,
                "tol"=>1e-5)
    TTCal.run(args)
    gains = npzread(gaintable)
    @test vecnorm(gains-truegains)/vecnorm(truegains) < 1e-5
end
test_one()

