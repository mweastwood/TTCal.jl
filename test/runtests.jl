using TTCal
using Base.Test
using Base.Dates
using CasaCore

function create_ms(name,x,y,z,ν)
    Nant  = length(x)
    Nfreq = length(ν)
    Nbase = div(Nant*(Nant-1),2) + Nant

    table = Table(name)
    addScalarColumn!(table,"ANTENNA1",Int32)
    addScalarColumn!(table,"ANTENNA2",Int32)
    addArrayColumn!(table,"UVW",Float64,[3])
    addScalarColumn!(table,"TIME",Float64)
    addArrayColumn!(table,"DATA",Complex64,[4,Nfreq])
    addArrayColumn!(table,"MODEL_DATA",Complex64,[4,Nfreq])
    addArrayColumn!(table,"CORRECTED_DATA",Complex64,[4,Nfreq])
    addRows!(table,Nbase)

    subtable = Table("$name/SPECTRAL_WINDOW")
    addArrayColumn!(subtable,"CHAN_FREQ",Cdouble,[Nfreq])
    addRows!(subtable,1)
    finalize(subtable)
    putKeyword!(table,"SPECTRAL_WINDOW","Table: $name/SPECTRAL_WINDOW")
    finalize(table)

    ms = MeasurementSet(name)
    # Initialize the Measurement Set with the interferometer details
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

    # Get a rough, current Julian date
    t = (2015.-1858.)*365.*24.*60.*60.

    putAntenna1!(ms,ant1)
    putAntenna2!(ms,ant2)
    putUVW!(ms,u,v,w)
    putTime!(ms,fill(t,Nbase))
    putFreq!(ms,ν)

    ms
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
    interferometer = TTCal.Interferometer(length(x),length(ν),1,zeros(Bool,0),zeros(Float64,length(x),2),0,0.)
    cyga = TTCal.Source("Cyg A",q"19h59m17.24s",q"+40d44m23.35s",21850.24,47e6,[-0.51,-0.18])
    sources = TTCal.getsources(ms)
    data = TTCal.visibilities(interferometer,ms,sources)
    putData!(ms,data)

    gains = TTCal.bandpass(interferometer,[ms],sources,TTCal.BandpassOptions(30,1e-5,4,true))
    truegains = ones(size(gains))

    # TODO: tighten this tolerance (at the moment I'm setting
    # it based on the tolerance passed to BandpassOptions)
    @test vecnorm(gains-truegains)/vecnorm(truegains) < 1e-5
end
test_one()

