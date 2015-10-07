using TTCal
using Base.Test
using CasaCore.Measures
using CasaCore.Tables

import TTCal.JonesMatrix

function xyz2uvw(x,y,z)
    Nant = length(x)
    Nbase = div(Nant*(Nant-1),2) + Nant
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

function ant1ant2(Nant)
    Nbase = div(Nant*(Nant-1),2) + Nant
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

function createms(Nant,Nfreq)
    Nbase = div(Nant*(Nant-1),2) + Nant

    x = 100*randn(Nant)
    y = 100*randn(Nant)
    z = randn(Nant)
    u,v,w = xyz2uvw(x,y,z)
    ν = linspace(40e6,60e6,Nfreq) |> collect
    t = (2015.-1858.)*365.*24.*60.*60. # a rough, current Julian date (in seconds)
    ant1,ant2 = ant1ant2(Nant)

    frame = ReferenceFrame()
    pos = observatory("OVRO_MMA")
    set!(frame,Epoch(epoch"UTC",Quantity(t,"s")))
    set!(frame,pos)
    zenith = Direction(dir"AZEL",q"0.0deg",q"90.0deg")
    phase_dir = measure(frame,zenith,dir"J2000")

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
    table["FLAG_ROW"] = zeros(Bool,Nbase)
    table["FLAG"] = zeros(Bool,4,Nfreq,Nbase)
    sources = readsources("sources.json")
    table["DATA"] = genvis(frame,phase_dir,sources,u,v,w,ν)

    name,table
end

srand(123)
include("jones.jl")
include("sourcemodel.jl")
include("fringepattern.jl")
include("subsrc.jl")
include("calibration.jl")
include("ampcal.jl")
include("gaincal.jl")
include("polcal.jl")
include("utm.jl")

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

