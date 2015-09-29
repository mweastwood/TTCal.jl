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

"""
These sources have a multi-component power-law spectrum
such that:

    log(flux) = log(I) + index[1]*log(ν/reffreq)
                       + index[2]*log²(ν/reffreq) + ...

Polarized fluxes are obtained in a similar manner by
substituting Q/U/V for I in the above expression.
"""
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
    @eval function $func{T<:AbstractFloat}(source::PointSource,frequency::T)
        powerlaw(source.$param,source.index,source.reffreq,frequency)
    end
    @eval function $func{T<:AbstractFloat}(source::PointSource,frequencies::Vector{T})
        [$func(source,frequency) for frequency in frequencies]
    end
end

flux{T<:AbstractFloat}(source::PointSource,frequency::T) = stokesI(source,frequency)
flux{T<:AbstractFloat}(source::PointSource,frequencies::Vector{T}) = stokesI(source,frequencies)

################################################################################
# Position

direction(source::PointSource) = source.dir

"""
Convert the direction into a right ascension and declination.
"""
function dir2radec(frame::ReferenceFrame,dir::Direction)
    j2000 = measure(frame,dir,Measures.J2000)
    ra  = longitude(j2000)
    dec = latitude(j2000)
    ra,dec
end

"""
Convert the direction into an azimuth and elevation.
"""
function dir2azel(frame::ReferenceFrame,dir::Direction)
    azel = measure(frame,dir,Measures.AZEL)
    az = longitude(azel)
    el = latitude(azel)
    az,el
end

"""
Convert the direction into the standard radio coordinate system.
"""
function dir2lm{ref}(frame::ReferenceFrame,phase_dir::Direction{Measures.J2000},dir::Direction{ref})
    long = longitude(phase_dir)
    lat  = latitude(phase_dir)
    θ1 = π/2 - long
    θ2 = π/2 - lat
    sin_θ1 = sin(θ1); cos_θ1 = cos(θ1)
    sin_θ2 = sin(θ2); cos_θ2 = cos(θ2)
    x,y,z = Measures.xyz_in_meters(measure(frame,dir,Measures.J2000))
    # - Rotate first by θ1 about the z-axis
    # - Then rotate by θ2 about the x-axis
    # ⌈    -l    ⌉   ⌈1    0        0    ⌉   ⌈+cos(θ1) -sin(θ1) 0⌉   ⌈x⌉
    # |    -m    | = |0 +cos(θ2) -sin(θ2)| * |+sin(θ1) +cos(θ1) 0| * |y|
    # ⌊√(1-l²-m²)⌋   ⌊0 +sin(θ2) +cos(θ2)⌋   ⌊   0        0     1⌋   ⌊z⌋
    l =        -cos_θ1*x +        sin_θ1*y
    m = -sin_θ1*cos_θ2*x - cos_θ1*cos_θ2*y + sin_θ2*z
    l,m
end

function lm2dir(phase_dir::Direction{Measures.J2000},l,m)
    long = longitude(phase_dir)
    lat  = latitude(phase_dir)
    θ1 = π/2 - long
    θ2 = π/2 - lat
    sin_θ1 = sin(θ1); cos_θ1 = cos(θ1)
    sin_θ2 = sin(θ2); cos_θ2 = cos(θ2)
    n = sqrt(1-l^2-m^2)
    # The rotation matrices in the previous function definition are easily
    # inverted by taking the transpose.
    # ⌈x⌉    ⌈+cos(θ1) +sin(θ1) 0⌉   ⌈1    0        0    ⌉   ⌈    -l    ⌉
    # |y| =  |-sin(θ1) +cos(θ1) 0| * |0 +cos(θ2) +sin(θ2)| * |    -m    |
    # ⌊z⌋    ⌊   0        0     1⌋   ⌊0 -sin(θ2) +cos(θ2)⌋   ⌊√(1-l²-m²)⌋
    x = -cos_θ1*l - sin_θ1*cos_θ2*m + sin_θ1*sin_θ2*n
    y = +sin_θ1*l - cos_θ1*cos_θ2*m + cos_θ1*sin_θ2*n
    z =                    sin_θ2*m +        cos_θ2*n
    Measures.from_xyz_in_meters(Measures.J2000,x,y,z)
end

radec(frame::ReferenceFrame,source::PointSource) = dir2radec(frame,direction(source))
azel(frame::ReferenceFrame,source::PointSource) = dir2azel(frame,direction(source))
lm(frame::ReferenceFrame,phase_dir::Direction,source::PointSource) = dir2lm(frame,phase_dir,direction(source))

function isabovehorizon(frame::ReferenceFrame,direction::Direction)
    az,el = dir2azel(frame,direction)
    el > 0
end

function isabovehorizon(frame::ReferenceFrame,source::PointSource)
    isabovehorizon(frame,direction(source))
end

################################################################################
# I/O

function readsources(filename::AbstractString)
    sources = PointSource[]
    parsed_sources = JSON.parsefile(filename)
    for parsed_source in parsed_sources
        name  = parsed_source["name"]
        if name == "Sun"
            dir = Direction(Measures.SUN)
        else
            ra  = parsed_source["ra"]
            dec = parsed_source["dec"]
            dir = Direction(Measures.J2000,Quanta.parse_ra(ra),Quanta.parse_dec(dec))
        end
        if haskey(parsed_source,"flux")
            I = parsed_source["flux"]
            Q = 0.0
            U = 0.0
            V = 0.0
        else
            I = parsed_source["I"]
            Q = parsed_source["Q"]
            U = parsed_source["U"]
            V = parsed_source["V"]
        end
        freq  = parsed_source["freq"]
        index = parsed_source["index"]
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
        if name != "Sun"
            dict["ra"]    = Quanta.format_ra(ra)
            dict["dec"]   = Quanta.format_dec(dec)
        end
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

