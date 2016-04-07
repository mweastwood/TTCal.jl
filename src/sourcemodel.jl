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

abstract Spectrum

function call(spectrum::Spectrum, ν::AbstractVector)
    StokesVector[spectrum(ν′) for ν′ in ν]
end

doc"""
    PowerLaw <: Spectrum

A multi-component power-law spectrum.

\\[
    \log_{10} S = \log_{10} S_0 + \sum_{n=1}^N \alpha_n \log_{10}\left(\frac{\nu}{\nu_0}\right)^n
\\]

where $S$ is a Stokes parameter, and $\alpha_n$ is the list of
spectral indices. At least one spectral index needs to be provided.
"""
type PowerLaw <: Spectrum
    stokes :: StokesVector
    ν :: Float64 # Hz
    α :: Vector{Float64}
end

PowerLaw(I,Q,U,V,ν,index) = PowerLaw(StokesVector(I,Q,U,V),ν,index)

function call(spectrum::PowerLaw, ν::AbstractFloat)
    s = sign(spectrum.stokes.I)
    log_I = log10(abs(spectrum.stokes.I))
    log_ν = log10(ν/spectrum.ν)
    for (i,α) in enumerate(spectrum.α)
        log_I += α*log_ν^i
    end
    I = s*10^log_I
    Q = spectrum.stokes.Q / spectrum.stokes.I * I
    U = spectrum.stokes.U / spectrum.stokes.I * I
    V = spectrum.stokes.V / spectrum.stokes.I * I
    StokesVector(I,Q,U,V)
end

function ==(lhs::PowerLaw, rhs::PowerLaw)
    lhs.stokes == rhs.stokes && lhs.ν == rhs.ν && lhs.α == rhs.α
end

type RFISpectrum <: Spectrum
    channels :: Vector{Float64}
    stokes   :: Vector{StokesVector}
end

function call(spectrum::RFISpectrum, ν::AbstractFloat)
    idx = searchsortedlast(spectrum.channels, ν)
    spectrum.stokes[idx]
end

abstract Source

"""
    PointSource <: Source

An astronomical point source.
"""
type PointSource <: Source
    name :: ASCIIString
    direction :: Direction
    spectrum  :: PowerLaw
end

"""
    GaussianSource <: Source

An astronomical Gaussian source.
"""
type GaussianSource <: Source
    name :: ASCIIString
    direction :: Direction
    spectrum  :: PowerLaw
    major_fwhm :: Float64 # FWHM along the major axis (radians)
    minor_fwhm :: Float64 # FWHM along the minor axis (radians)
    position_angle :: Float64 # (radians)
end

"""
    MultiSource <: Source

An astronomical source that has multiple components.
"""
type MultiSource <: Source
    name :: ASCIIString
    components :: Vector{Source}
end

"""
    RFISource <: Source

A terrestrial source of RFI. These sources are assumed to be
spectrally unsmooth and in the near field of the interferometer.
"""
type RFISource <: Source
    name :: ASCIIString
    position :: Position
    spectrum :: RFISpectrum
end

function ==(lhs::PointSource, rhs::PointSource)
    lhs.name == rhs.name && lhs.direction == rhs.direction && lhs.spectrum == rhs.spectrum
end

function ==(lhs::MultiSource, rhs::MultiSource)
    lhs.name == rhs.name && lhs.components == rhs.components
end

function isabovehorizon(frame::ReferenceFrame, direction::Direction)
    azel = measure(frame, direction, dir"AZEL")
    el = latitude(azel)
    el > 0
end

function isabovehorizon(frame::ReferenceFrame, source)
    isabovehorizon(frame, source.direction)
end

function isabovehorizon(frame::ReferenceFrame, source::MultiSource)
    for component in source.components
        if !isabovehorizon(frame, component)
            return false
        end
    end
    true
end

function abovehorizon{T<:Source}(frame::ReferenceFrame, sources::Vector{T})
    filter(sources) do source
        isabovehorizon(frame, source)
    end
end

"""
    readsources(filename)

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
    for c in parsed_sources
        source = construct_source(c)
        push!(sources, source)
    end
    sources
end

function construct_source(c)
    name = get(c, "name", "")
    if haskey(c, "components")
        # MultiSource
        components = Source[construct_source(dict) for dict in c["components"]]
        source = MultiSource(name, components)
    else
        dir  = get_source_direction(c)
        spec = get_source_spectrum(c)
        if haskey(c, "major-fwhm") && haskey(c, "minor-fwhm") && haskey(c, "position-angle")
            # GaussianSource
            major_fwhm = deg2rad(c["major-fwhm"]/3600)
            minor_fwhm = deg2rad(c["minor-fwhm"]/3600)
            position_angle = deg2rad(c["position-angle"])
            source = GaussianSource(name, dir, spec, major_fwhm, minor_fwhm, position_angle)
        else
            # PointSource
            source = PointSource(name, dir, spec)
        end
    end
    source
end

function get_source_direction(c)
    name = get(c, "name", "")
    if name == "Sun"
        dir = Direction(dir"SUN")
    elseif name == "Moon"
        dir = Direction(dir"MOON")
    elseif name == "Jupiter"
        dir = Direction(dir"JUPITER")
    else
        dir = Direction(dir"J2000", c["ra"], c["dec"])
    end
    dir
end

function get_source_spectrum(c)
    I = c["I"]
    Q = get(c, "Q", 0.0)
    U = get(c, "U", 0.0)
    V = get(c, "V", 0.0)
    freq  = c["freq"]
    index = c["index"]
    PowerLaw(I, Q, U, V, freq, index)
end

"""
    writesources(filename, sources)

Write the list of sources to the given location as a JSON file.
These sources can be read back in again using the `readsources` function.
"""
function writesources{T<:Source}(filename, sources::Vector{T})
    dicts = Dict{UTF8String,Any}[]
    for source in sources
        source_dict = deconstruct_source(source)
        push!(dicts, source_dict)
    end
    file = open(filename, "w")
    JSON.print(file, dicts)
    close(file)
    sources
end

function deconstruct_source(source::PointSource)
    c = Dict{UTF8String,Any}()
    c["name"] = source.name
    put_source_direction(c, source)
    put_source_spectrum(c, source)
    c
end

function deconstruct_source(source::GaussianSource)
    c = Dict{UTF8String,Any}()
    c["name"] = source.name
    put_source_direction(c, source)
    put_source_spectrum(c, source)
    c["major-fwhm"] = 3600*rad2deg(source.major_fwhm)
    c["minor-fwhm"] = 3600*rad2deg(source.minor_fwhm)
    c["position-angle"] = rad2deg(source.position-angle)
    c
end

function deconstruct_source(source::MultiSource)
    c = Dict{UTF8String,Any}()
    c["name"] = source.name
    c["components"] = [deconstruct_source(s) for s in source.components]
    c
end

function put_source_direction(c, source)
    if source.name != "Sun" && source.name != "Moon" && source.name != "Jupiter"
        ra  = longitude(source.direction) * radians
        dec =  latitude(source.direction) * radians
        c["ra"]  = sexagesimal(ra, hours = true)
        c["dec"] = sexagesimal(dec)
    end
end

function put_source_spectrum(c, source)
    c["I"]     = source.spectrum.stokes.I
    c["Q"]     = source.spectrum.stokes.Q
    c["U"]     = source.spectrum.stokes.U
    c["V"]     = source.spectrum.stokes.V
    c["freq"]  = source.spectrum.ν
    c["index"] = source.spectrum.α
end

