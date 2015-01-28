################################################################################
# Source Definitions

abstract AbstractSource

@doc """
These sources have a multi-component power-law spectrum
such that:

    log(S) = log(flux) + index[1]*log(ν/reffreq)
                       + index[2]*log²(ν/reffreq) + ...
""" ->
abstract PowerLawSource <: AbstractSource

@doc """
These are generic, run-of-the-mill sources that receive
a J2000 position and a power-law spectrum.
""" ->
type Source <: PowerLawSource
    name::ASCIIString
    dir::Direction
    flux::Float64
    reffreq::quantity(Float64,Hertz)
    index::Vector{Float64}
end

@doc """
Solar system objects move through the sky and hence their
J2000 position must be calculated for the given epoch.
""" ->
type SolarSystemSource <: PowerLawSource
    name::ASCIIString
    flux::Float64
    reffreq::quantity(Float64,Hertz)
    index::Vector{Float64}
end

@doc """
RFI usually has an unsmooth spectrum and does not rotate
with the sky. Hence these sources are allowed to have any
arbitrary spectrum, and receive an AZEL position.
""" ->
type RFISource <: AbstractSource
    name::ASCIIString
    dir::Direction
    flux::Vector{Float64}
    freq::Vector{quantity(Float64,Hertz)}
end

function RFISource(frame::ReferenceFrame,
                   source::Source,
                   frequencies::Vector{quantity(Float64,Hertz)})
    dir = measure(frame,source.dir,"AZEL")
    S = flux(source,frequencies)
    RFISource(source.name,dir,S,frequencies)
end

################################################################################
# Flux

function flux(source::PowerLawSource,frequency::quantity(Float64,Hertz))
    logflux = log10(source.flux)
    logfreq = log10(frequency/source.reffreq)
    for (i,index) in enumerate(source.index)
        logflux += index*logfreq.^i
    end
    10.0.^logflux
end

function flux(source::RFISource,frequencies::Vector{quantity(Float64,Hertz)})
    frequencies == source.freq || error("Provided list of frequencies does not match.")
    source.flux
end

flux(source::RFISource) = source.flux

function flux(source::PowerLawSource,frequencies::Vector{quantity(Float64,Hertz)})
    [flux(source,frequency) for frequency in frequencies]
end

################################################################################
# Position

direction(source::Source) = source.dir
direction(source::SolarSystemSource) = Direction(source.name)
direction(source::RFISource) = source.dir

@doc """
Convert the direction into an azimuth and elevation.
""" ->
function dir2azel(frame::ReferenceFrame,dir::Direction)
    dir = measure(frame,dir,"AZEL")
    az = dir.m[1]
    el = dir.m[2]
    az,el
end

@doc """
Convert the direction into the standard radio coordinate system.
""" ->
function dir2lm(frame::ReferenceFrame,dir::Direction)
    az,el = dir2azel(frame,dir)
    azel2lm(az,el)
end

function azel2lm(az,el)
    l = cos(el)*sin(az)
    m = cos(el)*cos(az)
    l,m
end

function lm2azel(l,m)
    az = atan2(l,m)*Radian
    el = acos(sqrt(l^2+m^2))*Radian
    az,el
end

azel(frame::ReferenceFrame,source::AbstractSource) = dir2azel(frame,direction(source))
lm(frame::ReferenceFrame,source::AbstractSource) = dir2lm(frame,direction(source))

@doc """
Returns true if the source is above the horizon, false if the source
is below the horizon.
""" ->
function isabovehorizon(frame::ReferenceFrame,source::Source)
    az,el = azel(frame,source)
    ifelse(el > 0.0Radian,true,false)
end

################################################################################
# I/O

function readsources(filename::AbstractString)
    sources = TTCal.Source[]
    parsed_sources = JSON.parsefile(filename)
    for parsed_source in parsed_sources
        name  = parsed_source["name"]
        ra    = parsed_source["ra"]
        dec   = parsed_source["dec"]
        flux  = parsed_source["flux"]
        freq  = parsed_source["freq"]
        index = parsed_source["index"]
        dir = Direction("J2000",ra_str(ra),dec_str(dec))
        push!(sources,TTCal.Source(name,dir,flux,freq*Hertz,index))
    end
    sources
end

function writesources(filename::AbstractString,sources::Vector{Source})
    dicts = Dict{UTF8String,Any}[]
    for source in sources
        dict = Dict{UTF8String,Any}()
        dict["ref"]   = "TTCal"
        dict["name"]  = source.name
        dict["ra"]    =  ra_str(source.dir.m[1])
        dict["dec"]   = dec_str(source.dir.m[2])
        dict["flux"]  = source.flux
        dict["freq"]  = float(source.reffreq)
        dict["index"] = source.index
        push!(dicts,dict)
    end
    file = open(filename,"w")
    JSON.print(file,dicts)
    close(file)
    nothing
end

