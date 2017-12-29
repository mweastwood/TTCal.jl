# Copyright (c) 2015-2017 Michael Eastwood
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

struct SkyModel
    sources :: Vector{Source}
end

function SkyModel(source::Source)
    SkyModel([source])
end

"""
    readsky(filename)

Read the list of point sources from the given JSON file.  The format must be as follows:

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

A right ascension and declination does not need to be specified if the name of the source is "Sun",
"Moon", or "Jupiter". These sources will have their location automatically determined by CasaCore.
"""
function readsky(filename)
    sources = Source[]
    parsed_sources = JSON.parsefile(filename)
    for c in parsed_sources
        source = construct_source(c)
        push!(sources, source)
    end
    SkyModel(sources)
end

function construct_source(c)
    name   = get(c, "name", "")
    shapes = get_shapes(name, c)
    Source(name, shapes)
end

function get_shapes(name, c)
    if haskey(c, "components")
        return [get_shapes(name, d) for d in c["components"]]
    #elseif haskey(c, "rfi-frequencies") && haskey(c, "rfi-I")
    #    # RFISource
    #    position = get_source_position(c)
    #    spectrum = get_rfi_spectrum(c)
    #    source = RFISource(name, position, spectrum)
    else
        dir  = get_source_direction(name, c)
        spec = get_source_spectrum(c)
        if haskey(c, "major-fwhm") && haskey(c, "minor-fwhm") && haskey(c, "position-angle")
            # Gaussian
            major_fwhm = deg2rad(c["major-fwhm"]/3600)
            minor_fwhm = deg2rad(c["minor-fwhm"]/3600)
            position_angle = deg2rad(c["position-angle"])
            return Gaussian(dir, spec, major_fwhm, minor_fwhm, position_angle)
        #elseif haskey(c, "radius")
        #    # Disk
        #    radius = deg2rad(c["radius"]/3600)
        #    source = DiskSource(name, dir, spec, radius)
        elseif haskey(c, "coefficients")
            # Shapelet
            scale = c["scale-parameter"]
            coeff = convert(Vector{Float64}, c["coefficients"])
            source = Shapelet(dir, spec, scale, coeff)
        else
            # Point
            return Point(dir, spec)
        end
    end
end

function get_source_direction(name, c)
    if name == "Sun"
        dir = Direction(dir"SUN")
    elseif name == "Moon"
        dir = Direction(dir"MOON")
    elseif name == "Jupiter"
        dir = Direction(dir"JUPITER")
    elseif haskey(c, "az") && haskey(c, "el")
        dir = Direction(dir"AZEL", c["az"], c["el"])
    else
        dir = Direction(dir"J2000", c["ra"], c["dec"])
    end
    dir
end

#function get_source_position(c)
#    sys  = c["sys"]
#    el   = c["el"]*meters
#    long = c["long"]*degrees
#    lat  = c["lat"]*degrees
#    if sys == "WGS84"
#        pos = Position(pos"WGS84", el, long, lat)
#    elseif sys == "ITRF"
#        pos = Position(pos"ITRF", el, long, lat)
#    else
#        error("unknown coordinate system")
#    end
#    pos
#end

function get_source_spectrum(c)
    I = c["I"]
    Q = get(c, "Q", 0.0)
    U = get(c, "U", 0.0)
    V = get(c, "V", 0.0)
    freq  = c["freq"]
    index = c["index"]
    PowerLaw(I, Q, U, V, freq*u"Hz", index)
end

#function get_rfi_spectrum(c)
#    channels = c["rfi-frequencies"]
#    N = length(channels)
#    I = c["rfi-I"]
#    Q = get(c, "rfi-Q", zeros(N))
#    U = get(c, "rfi-U", zeros(N))
#    V = get(c, "rfi-V", zeros(N))
#    stokes = [StokesVector(I[idx], Q[idx], U[idx], V[idx]) for idx = 1:N]
#    RFISpectrum(channels, stokes)
#end
#
#"""
#    writesources(filename, sources)
#
#Write the list of sources to the given location as a JSON file.  These sources can be read back in
#again using the `readsky` function.
#"""
#function writesources{T<:Source}(filename, sources::Vector{T})
#    dicts = Dict{String,Any}[]
#    for source in sources
#        source_dict = deconstruct_source(source)
#        push!(dicts, source_dict)
#    end
#    file = open(filename, "w")
#    JSON.print(file, dicts)
#    close(file)
#    sources
#end
#
#function deconstruct_source(source::PointSource)
#    c = Dict{String,Any}()
#    c["name"] = source.name
#    put_source_direction(c, source)
#    put_source_spectrum(c, source)
#    c
#end
#
#function deconstruct_source(source::GaussianSource)
#    c = Dict{String,Any}()
#    c["name"] = source.name
#    put_source_direction(c, source)
#    put_source_spectrum(c, source)
#    c["major-fwhm"] = 3600*rad2deg(source.major_fwhm)
#    c["minor-fwhm"] = 3600*rad2deg(source.minor_fwhm)
#    c["position-angle"] = rad2deg(source.position_angle)
#    c
#end
#
#function deconstruct_source(source::ShapeletSource)
#    c = Dict{String,Any}()
#    c["name"] = source.name
#    put_source_direction(c, source)
#    put_source_spectrum(c, source)
#    c["scale-parameter"] = source.scale
#    c["coefficients"] = source.coeff
#    c
#end
#
#function deconstruct_source(source::MultiSource)
#    c = Dict{String,Any}()
#    c["name"] = source.name
#    c["components"] = [deconstruct_source(s) for s in source.components]
#    c
#end
#
#function deconstruct_source(source::RFISource)
#    c = Dict{String,Any}()
#    c["name"] = source.name
#    put_source_position(c, source)
#    put_rfi_spectrum(c, source)
#    c
#end
#
#function put_source_direction(c, source)
#    if source.name != "Sun" && source.name != "Moon" && source.name != "Jupiter"
#        ra  = longitude(source.direction) * radians
#        dec =  latitude(source.direction) * radians
#        c["ra"]  = sexagesimal(ra, digits = 3, hours = true)
#        c["dec"] = sexagesimal(dec, digits = 2)
#    end
#end
#
#function put_source_position(c, source)
#    pos = source.position
#    if pos.sys === pos"WGS84"
#        sys = "WGS84"
#    elseif pos.sys === pos"ITRF"
#        sys = "ITRF"
#    else
#        error("unknown coordinate system")
#    end
#    el   =    radius(pos)
#    long = longitude(pos) |> rad2deg
#    lat  =  latitude(pos) |> rad2deg
#    c["sys"]  = sys
#    c["el"]   = el
#    c["long"] = long
#    c["lat"]  = lat
#end
#
#function put_source_spectrum(c, source)
#    c["I"]     = source.spectrum.stokes.I
#    c["Q"]     = source.spectrum.stokes.Q
#    c["U"]     = source.spectrum.stokes.U
#    c["V"]     = source.spectrum.stokes.V
#    c["freq"]  = source.spectrum.ν
#    c["index"] = source.spectrum.α
#end
#
#function put_rfi_spectrum(c, source)
#    spec = source.spectrum
#    channels = spec.channels
#    stokes = spec.stokes
#    N = length(channels)
#    I = Float64[stokes[idx].I for idx = 1:N]
#    Q = Float64[stokes[idx].Q for idx = 1:N]
#    U = Float64[stokes[idx].U for idx = 1:N]
#    V = Float64[stokes[idx].V for idx = 1:N]
#    c["rfi-frequencies"] = channels
#    c["rfi-I"] = I
#    c["rfi-Q"] = Q
#    c["rfi-U"] = U
#    c["rfi-V"] = V
#end

