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
    reffreq::Float64
    index::Vector{Float64}
end

name(source::PointSource) = source.name
Base.show(io::IO,source::PointSource) = print(io,name(source))

################################################################################
# Flux

function powerlaw(reference_flux,index,reference_frequency,frequency)
    s = sign(reference_flux)
    logflux = log10(abs(reference_flux))
    logfreq = log10(frequency/reference_frequency)
    for (i,α) in enumerate(index)
        logflux += α*logfreq.^i
    end
    s*10.0.^logflux
end

for param in (:I,:Q,:U,:V)
    func = symbol("stokes",param)
    @eval function $func{T<:FloatingPoint}(source::PointSource,frequency::T)
        powerlaw(source.$param,source.index,source.reffreq,frequency)
    end
    @eval function $func{T<:FloatingPoint}(source::PointSource,frequencies::Vector{T})
        [$func(source,frequency) for frequency in frequencies]
    end
end

flux{T<:FloatingPoint}(source::PointSource,frequency::T) = stokesI(source,frequency)
flux{T<:FloatingPoint}(source::PointSource,frequencies::Vector{T}) = stokesI(source,frequencies)

################################################################################
# Position

direction(source::PointSource) = source.dir

@doc """
Convert the direction into a right ascension and declination.
""" ->
function dir2radec(frame::ReferenceFrame,dir::Direction)
    j2000 = measure(frame,dir,Measures.J2000)
    ra  = longitude(azel)
    dec = latitude(azel)
    ra,dec
end

@doc """
Convert the direction into an azimuth and elevation.
""" ->
function dir2azel(frame::ReferenceFrame,dir::Direction)
    azel = measure(frame,dir,Measures.AZEL)
    az = longitude(azel)
    el = latitude(azel)
    az,el
end

@doc """
Convert the direction into the standard radio coordinate system.
""" ->
function dir2lm{ref}(phase_dir::Direction{ref},dir::Direction{ref})
    long = longitude(phase_dir)
    lat  = latitude(phase_dir)
    θ1 = π/2 - long
    θ2 = π/2 - lat
    sin_θ1 = sin(θ1); cos_θ1 = cos(θ1)
    sin_θ2 = sin(θ2); cos_θ2 = cos(θ2)
    x,y,z = Measures.xyz_in_meters(dir)
    # - Rotate first by θ1 about the z-axis
    # - Then rotate by θ2 about the x-axis
    # ⌈    -l    ⌉   ⌈1    0        0    ⌉   ⌈+cos(θ1) -sin(θ1) 0⌉   ⌈x⌉
    # |    -m    | = |0 +cos(θ2) -sin(θ2)| * |+sin(θ1) +cos(θ1) 0| * |y|
    # ⌊√(1-l²-m²)⌋   ⌊0 +sin(θ2) +cos(θ2)⌋   ⌊   0        0     1⌋   ⌊z⌋
    l =       -cos_θ1*x  +        sin_θ1*y
    m = -sin_θ1*cos_θ2*x - cos_θ1*cos_θ2*y + sin_θ2*z
    l,m
end

function azel2lm(az,el)
    l = cos(el).*sin(az)
    m = cos(el).*cos(az)
    l,m
end

function lm2azel(l,m)
    az = atan2(l,m)
    el = acos(hypot(l,m))
    az,el
end

radec(frame::ReferenceFrame,source::PointSource) = dir2radec(frame,direction(source))
azel(frame::ReferenceFrame,source::PointSource) = dir2azel(frame,direction(source))
lm(phase_dir::Direction,source::PointSource) = dir2lm(phase_dir,direction(source))

function isabovehorizon(frame::ReferenceFrame,source::PointSource)
    az,el = azel(frame,source)
    el > 0
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
        freq  = parsed_source["freq"]
        index = parsed_source["index"]
        dir = Direction(Measures.J2000,Quanta.parse_ra(ra),Quanta.parse_dec(dec))
        push!(sources,PointSource(name,dir,I,Q,U,V,freq,index))
    end
    sources
end

function writesources(filename::AbstractString,sources::Vector{PointSource})
    dicts = Dict{UTF8String,Any}[]
    for source in sources
        ra  = longitude(source.dir,Degree)
        dec =  latitude(source.dir,Degree)
        dict = Dict{UTF8String,Any}()
        dict["ref"]   = "TTCal"
        dict["name"]  = source.name
        dict["ra"]    = Quanta.format_ra(ra)
        dict["dec"]   = Quanta.format_dec(dec)
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

