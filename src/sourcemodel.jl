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

"""
    type PointSource

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
direction(source::PointSource) = source.dir

function powerlaw(reference_flux,index,reference_frequency,frequency)
    s = sign(reference_flux)
    logflux = log10(abs(reference_flux))
    logfreq = log10(frequency/reference_frequency)
    for (i,α) in enumerate(index)
        logflux += α*logfreq.^i
    end
    s*10.0.^logflux
end

function flux(source::PointSource, frequency::Number)
    powerlaw(source.I,source.index,source.reffreq,frequency)
end

function flux{T<:Number}(source::PointSource, frequencies::Vector{T})
    T[flux(source,frequency) for frequency in frequencies]
end

function stokes_flux(source::PointSource, frequency::Number)
    I = powerlaw(source.I,source.index,source.reffreq,frequency)
    Q = powerlaw(source.Q,source.index,source.reffreq,frequency)
    U = powerlaw(source.U,source.index,source.reffreq,frequency)
    V = powerlaw(source.V,source.index,source.reffreq,frequency)
    [I,Q,U,V]
end

"""
Convert the direction into a right ascension and declination.
"""
function dir2radec(frame::ReferenceFrame,dir::Direction)
    j2000 = measure(frame,dir,dir"J2000")
    ra  = longitude(j2000)
    dec =  latitude(j2000)
    ra,dec
end

"""
Convert the direction into an azimuth and elevation.
"""
function dir2azel(frame::ReferenceFrame,dir::Direction)
    azel = measure(frame,dir,dir"AZEL")
    az = longitude(azel)
    el =  latitude(azel)
    az,el
end

"""
Convert the direction into the standard radio coordinate system.
"""
function dir2lm{ref}(frame::ReferenceFrame,phase_dir::Direction{dir"J2000"},dir::Direction{ref})
    long = longitude(phase_dir)
    lat  = latitude(phase_dir)
    θ1 = π/2 - long
    θ2 = π/2 - lat
    sin_θ1 = sin(θ1); cos_θ1 = cos(θ1)
    sin_θ2 = sin(θ2); cos_θ2 = cos(θ2)
    x,y,z = Measures.xyz_in_meters(measure(frame,dir,dir"J2000"))
    # - Rotate first by θ1 about the z-axis
    # - Then rotate by θ2 about the x-axis
    # ⌈    -l    ⌉   ⌈1    0        0    ⌉   ⌈+cos(θ1) -sin(θ1) 0⌉   ⌈x⌉
    # |    -m    | = |0 +cos(θ2) -sin(θ2)| * |+sin(θ1) +cos(θ1) 0| * |y|
    # ⌊√(1-l²-m²)⌋   ⌊0 +sin(θ2) +cos(θ2)⌋   ⌊   0        0     1⌋   ⌊z⌋
    l =        -cos_θ1*x +        sin_θ1*y
    m = -sin_θ1*cos_θ2*x - cos_θ1*cos_θ2*y + sin_θ2*z
    l,m
end

function lm2dir(phase_dir::Direction{dir"J2000"},l,m)
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
    Measures.from_xyz_in_meters(dir"J2000",x,y,z)
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
            dir = Direction(dir"SUN")
        else
            ra  = parsed_source["ra"]
            dec = parsed_source["dec"]
            dir = Direction(dir"J2000",Quantity(Measures.sexagesimal(ra),"deg"),
                                       Quantity(Measures.sexagesimal(dec),"deg"))
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
        ra  = longitude(source.dir,"deg")
        dec =  latitude(source.dir,"deg")
        dict = Dict{UTF8String,Any}()
        dict["ref"]   = "TTCal"
        dict["name"]  = source.name
        if name != "Sun"
            dict["ra"]  = format_ra(ra)
            dict["dec"] = format_dec(dec)
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

"""
    write_ds9_regions(filename,sources::Vector{PointSource})

Write the list of sources to a DS9 region file. This file can then be loaded
into DS9 to highlight sources within images.
"""
function write_ds9_regions(filename,sources::Vector{PointSource})
    open(filename,"w") do f
        Base.write(f,"global color=red edit=0 move=0 delete=1\n")
        Base.write(f,"fk5\n")
        for source in sources
            source.name == "Sun" && continue
            ra  = longitude(source.dir,"deg")
            dec =  latitude(source.dir,"deg")
            Base.write(f,"circle($(format_ra(ra)),$(format_dec(dec)),1000\") # text = {$(source.name)}\n")
        end
    end
end

function format_ra(ra::Float64)
    ra /= 15
    ra  = mod(ra,24)
    hrs = floor(Integer,ra)
    ra  = (ra-hrs)*60
    min = floor(Integer,ra)
    ra  = (ra-min)*60
    sec = ra
    @sprintf("%dh%02dm%07.4fs",hrs,min,sec)
end

function format_dec(dec::Float64)
    s = sign(dec)
    dec *= s
    dec = mod(dec,90)
    deg = floor(Integer,dec)
    dec = (dec-deg)*60
    min = floor(Integer,dec)
    dec = (dec-min)*60
    sec = dec
    @sprintf("%+dd%02dm%07.4fs",s*deg,min,sec)
end

