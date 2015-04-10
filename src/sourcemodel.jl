# Copyright (c) 2015 Michael Eastwood
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

################################################################################
# Source Definitions

@doc """
These sources have a multi-component power-law spectrum
such that:

    log(flux) = log(I) + index[1]*log(ν/reffreq)
                       + index[2]*log²(ν/reffreq) + ...

Polarized fluxes are obtained in a similar manner by
substituting Q/U/V for I in the above expression.
""" ->
type PointSource
    name::ASCIIString
    dir::Direction
    I::Float64
    Q::Float64
    U::Float64
    V::Float64
    reffreq::quantity(Float64,Hertz)
    index::Vector{Float64}
end

Base.show(io::IO,source::PointSource) = print(io,source.name)

################################################################################
# Flux

function getflux(reference_flux,index,reference_frequency,frequency)
    s = sign(reference_flux)
    logflux = log10(abs(reference_flux))
    logfreq = log10(frequency/reference_frequency)
    for (i,α) in enumerate(index)
        logflux += α*logfreq.^i
    end
    s*10.0.^logflux
end

for param in (:I,:Q,:U,:V)
    name = symbol("getstokes",param)
    @eval function $name{T<:FloatingPoint}(source::PointSource,frequency::quantity(T,Hertz))
        getflux(source.$param,source.index,source.reffreq,frequency)
    end
    @eval function $name{T<:FloatingPoint}(source::PointSource,frequencies::Vector{quantity(T,Hertz)})
        [$name(source,frequency) for frequency in frequencies]
    end
end

getflux{T<:FloatingPoint}(source::PointSource,frequency::quantity(T,Hertz)) = getstokesI(source,frequency)
getflux{T<:FloatingPoint}(source::PointSource,frequencies::Vector{quantity(T,Hertz)}) = getstokesI(source,frequencies)

################################################################################
# Position

direction(source::PointSource) = source.dir

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
    l = cos(el).*sin(az)
    m = cos(el).*cos(az)
    l,m
end

function lm2azel(l,m)
    az = atan2(l,m)*Radian
    el = acos(sqrt(l^2+m^2))*Radian
    az,el
end

getazel(frame::ReferenceFrame,source::PointSource) = dir2azel(frame,direction(source))
getlm(frame::ReferenceFrame,source::PointSource) = dir2lm(frame,direction(source))

@doc """
Returns true if the source is above the horizon, false if the source
is below the horizon.
""" ->
isabovehorizon(frame::ReferenceFrame,source::PointSource) = isabovehorizon(frame,direction(source))

function isabovehorizon(frame::ReferenceFrame,dir::Direction)
    az,el = dir2azel(frame,dir)
    ifelse(el > 0.0Radian,true,false)
end

################################################################################
# I/O

function readsources(filename::AbstractString)
    sources = PointSource[]
    parsed_sources = JSON.parsefile(filename)
    for parsed_source in parsed_sources
        name  = parsed_source["name"]
        ra    = parsed_source["ra"]
        dec   = parsed_source["dec"]
        I     = parsed_source["I"]
        Q     = parsed_source["Q"]
        U     = parsed_source["U"]
        V     = parsed_source["V"]
        freq  = parsed_source["freq"]*Hertz
        index = parsed_source["index"]
        dir = Direction("J2000",ra_str(ra),dec_str(dec))
        push!(sources,PointSource(name,dir,I,Q,U,V,freq,index))
    end
    sources
end

function writesources(filename::AbstractString,sources::Vector{PointSource})
    dicts = Dict{UTF8String,Any}[]
    for source in sources
        dict = Dict{UTF8String,Any}()
        dict["ref"]   = "TTCal"
        dict["name"]  = source.name
        dict["ra"]    =  ra_str(source.dir.m[1])
        dict["dec"]   = dec_str(source.dir.m[2])
        dict["I"]     = source.I
        dict["Q"]     = source.Q
        dict["U"]     = source.U
        dict["V"]     = source.V
        dict["freq"]  = float(source.reffreq)
        dict["index"] = source.index
        push!(dicts,dict)
    end
    file = open(filename,"w")
    JSON.print(file,dicts)
    close(file)
    nothing
end

