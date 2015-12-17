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

doc"""
    immutable Spectrum

These sources have a multi-component power-law spectrum
such that:

\\[
    \log_{10} S = \log_{10} S_0 + \sum_{n=1}^N \alpha_n \log_{10}\left(\frac{\nu}{\nu_0}\right)^n
\\]

where $S$ is a Stokes parameter, and $\alpha_n$ is the list of
spectral indices. At least one spectral index needs to be provided.
"""
immutable Spectrum
    stokes::StokesVector
    ν0::Float64 # Hz
    spectral_index::Vector{Float64}
end

Spectrum(I,Q,U,V,ν,index) = Spectrum(StokesVector(I,Q,U,V),ν,index)

function call(spectrum::Spectrum, ν::AbstractFloat)
    s = sign(spectrum.stokes.I)
    log_I = log10(abs(spectrum.stokes.I))
    log_ν = log10(ν/spectrum.ν0)
    for (i,α) in enumerate(spectrum.spectral_index)
        log_I += α*log_ν^i
    end
    I = s*10^log_I
    Q = spectrum.stokes.Q / spectrum.stokes.I * I
    U = spectrum.stokes.U / spectrum.stokes.I * I
    V = spectrum.stokes.V / spectrum.stokes.I * I
    StokesVector(I,Q,U,V)
end

function call(spectrum::Spectrum, ν::AbstractVector)
    [spectrum(ν[β]) for β = 1:length(ν)]
end

abstract Component

"""
    immutable Source

This type represents a radio source.

Each source is composed of one or more components.
"""
immutable Source
    name::ASCIIString
    components::Vector{Component}
end

Source(name::ASCIIString,component::Component) = Source(name,[component])

immutable Point <: Component
    name::ASCIIString
    direction::Direction
    spectrum::Spectrum
end

# TODO
#immutable Gaussian <: Component
#end

# TODO
#immutable Ellipse <: Component
#end

"""
    j2000_radec(frame::ReferenceFrame, component::Component) -> ra,dec

Compute the J2000 right ascension and declination of the component (in radians).
"""
j2000_radec(frame::ReferenceFrame, component::Component) = j2000_radec(frame,component.direction)

function j2000_radec(frame::ReferenceFrame, direction::Direction)
    j2000 = measure(frame,direction,dir"J2000")
    ra  = longitude(j2000)
    dec =  latitude(j2000)
    ra,dec
end

"""
    local_azel(frame::ReferenceFrame, component::Component) -> az,el

Compute the local azimuth and elevation of the component (in radians).
"""
local_azel(frame::ReferenceFrame, component::Component) = local_azel(frame,component.direction)

function local_azel(frame::ReferenceFrame, direction::Direction)
    azel = measure(frame,direction,dir"AZEL")
    az = longitude(azel)
    el =  latitude(azel)
    az,el
end

function isabovehorizon(frame::ReferenceFrame, direction::Direction)
    az,el = local_azel(frame,direction)
    el > 0
end

function isabovehorizon(frame::ReferenceFrame, component::Point)
    isabovehorizon(frame,component.direction)
end

function isabovehorizon(frame::ReferenceFrame, source::Source)
    for component in source.components
        if !isabovehorizon(frame,component)
            return false
        end
    end
    true
end

function abovehorizon(frame::ReferenceFrame, sources::Vector{Source})
    filter(sources) do source
        isabovehorizon(frame,source)
    end
end

doc"""
    direction_cosines(phase_dir::Direction, dir::Direction) -> l,m

Compute the direction cosines $(l,m)$ for the given direction with respect to the
phase direction.
"""
function direction_cosines(phase_dir::Direction, dir::Direction)
    long = longitude(phase_dir)
    lat  =  latitude(phase_dir)
    θ1 = π/2 - long
    θ2 = π/2 - lat
    sin_θ1 = sin(θ1); cos_θ1 = cos(θ1)
    sin_θ2 = sin(θ2); cos_θ2 = cos(θ2)
    x,y,z = dir.x, dir.y, dir.z
    # Rotate first by θ1 about the z-axis
    # Then rotate by θ2 about the x-axis
    # ⌈    -l    ⌉   ⌈1    0        0    ⌉   ⌈+cos(θ1) -sin(θ1) 0⌉   ⌈x⌉
    # |    -m    | = |0 +cos(θ2) -sin(θ2)| * |+sin(θ1) +cos(θ1) 0| * |y|
    # ⌊√(1-l²-m²)⌋   ⌊0 +sin(θ2) +cos(θ2)⌋   ⌊   0        0     1⌋   ⌊z⌋
    l =        -cos_θ1*x +        sin_θ1*y
    m = -sin_θ1*cos_θ2*x - cos_θ1*cos_θ2*y + sin_θ2*z
    l,m
end

function undo_direction_cosines(phase_dir::Direction,l,m)
    long = longitude(phase_dir)
    lat  =  latitude(phase_dir)
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
    Direction(dir"J2000",x,y,z)
end

"""
    readsources(filename) -> Vector{Source}

Read the list of point sources from the given JSON file.
The format must be as follows:

    [
        {
            "ref": "Baars et al. 1977",
            "name": "Cas A",
            "ra": "23h23m24s",
            "dec": "58d48m54s",
            "I": 555904.26,
            "freq": 1.0e6,
            "index": [-0.770]
        },
        {
            "ref": "Baars et al. 1977",
            "name": "Cyg A",
            "ra": "19h59m28.35663s",
            "dec": "+40d44m02.0970s",
            "I": 49545.02,
            "freq": 1.0e6,
            "index": [+0.085,-0.178]
        }
    ]

Additional Stokes parameters may also be specified.

    {
        ...
        "I": 100
        "Q": 20
        "U": 3.14
        "V": -30
        ...
    }

A right ascension and declination does not need to be specified if
the name of the source is "Sun", "Moon", or "Jupiter". These sources
will have their location automatically determined by CasaCore.
"""
function readsources(filename)
    sources = Source[]
    parsed_sources = JSON.parsefile(filename)
    for parsed_source in parsed_sources
        name = parsed_source["name"]
        components = Component[]
        if haskey(parsed_source,"components")
            parsed_components = parsed_source["components"]
            for parsed_component in parsed_components
                component = construct_component(parsed_component)
                push!(components,component)
            end
        else
            component = construct_component(parsed_source)
            push!(components,component)
        end
        push!(sources,Source(name,components))
    end
    sources
end

function construct_component(c)
    name = get(c,"name","")

    if name == "Sun"
        dir = Direction(dir"SUN")
    elseif name == "Moon"
        dir = Direction(dir"MOON")
    elseif name == "Jupiter"
        dir = Direction(dir"JUPITER")
    else
        dir = Direction(dir"J2000", c["ra"], c["dec"])
    end

    if haskey(c,"flux")
        warn("""
            $filename is out of date
            Replace "flux" with "I" (for the Stokes I flux).
            Additional entries for "Q", "U", and "V" may also be added.
        """)
        # for compatibility with old sources.json files,
        # which used "flux" in place of the Stokes parameters.
        I = c["flux"]
    else
        I = c["I"]
    end
    Q = get(c,"Q",0.0)
    U = get(c,"U",0.0)
    V = get(c,"V",0.0)

    freq  = c["freq"]
    index = c["index"]

    Point(name,dir,Spectrum(I,Q,U,V,freq,index))
end

"""
    writesources(filename, sources::Vector{Source})

Write the list of sources to the given location as a JSON file.
These sources can be read back in again using the `readsources` function.
"""
function writesources(filename, sources::Vector{Source})
    dicts = Dict{UTF8String,Any}[]
    for source in sources
        source_dict = Dict{UTF8String,Any}()
        source_dict["ref"]  = "TTCal"
        source_dict["name"] = source.name

        component_dicts = Dict{UTF8String,Any}[]
        for (n,component) in enumerate(source.components)
            component_dict = Dict{UTF8String,Any}()
            component_dict["name"] = component.name

            if component.name != "Sun" && component.name != "Moon" && component.name != "Jupiter"
                ra  = longitude(component.direction) * radians
                dec =  latitude(component.direction) * radians
                component_dict["ra"]  = sexagesimal(ra, hours = true)
                component_dict["dec"] = sexagesimal(dec)
            end

            component_dict["I"]     = component.spectrum.stokes.I
            component_dict["Q"]     = component.spectrum.stokes.Q
            component_dict["U"]     = component.spectrum.stokes.U
            component_dict["V"]     = component.spectrum.stokes.V
            component_dict["freq"]  = component.spectrum.ν0
            component_dict["index"] = component.spectrum.spectral_index

            push!(component_dicts,component_dict)
        end

        source_dict["components"] = component_dicts
        push!(dicts,source_dict)
    end
    file = open(filename,"w")
    JSON.print(file,dicts)
    close(file)
    sources
end

#=
"""
    write_ds9_regions(filename, sources::Vector{Source})

Write the list of sources to a DS9 region file. This file can then be loaded
into DS9 to highlight sources within images.
"""
function write_ds9_regions(filename,sources::Vector{Source})
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
=#

