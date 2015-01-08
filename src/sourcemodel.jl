# The source flux is modeled as follows:
#     log(S) = log(flux) + index[1]*log(ν/reffreq)
#                        + index[2]*log²(ν/reffreq) + ...
type Source
    name::ASCIIString
    dir::Direction
    flux::Float64
    reffreq::quantity(Float64,Hertz)
    index::Vector{Float64}
end

Source(name,ra,dec,flux,reffreq,index::Float64) = Source(name,ra,dec,flux,reffreq,[index])

function getflux(source::Source,frequency::quantity(Float64,Hertz))
    logflux = log10(source.flux)
    for (i,index) in enumerate(source.index)
        logflux += index*log10(frequency/source.reffreq).^i
    end
    10.0.^logflux
end

function getflux(source::Source,frequencies::Vector{quantity(Float64,Hertz)})
    [getflux(source,frequency) for frequency in frequencies]
end

function getazel(frame::ReferenceFrame,dir::Direction)
    if dir.system != "AZEL"
        dir = measure(frame,dir,"AZEL")
    end
    az = dir.m[1]
    el = dir.m[2]
    az,el
end
getazel(frame::ReferenceFrame,ra,dec) = getazel(frame,Direction("J2000",ra,dec))
getazel(frame::ReferenceFrame,source::Source) = getazel(frame,source.dir)

@doc """
Convert a given RA and dec to the standard radio coordinate system.
""" ->
function getlm(frame::ReferenceFrame,dir::Direction)
    az,el = getazel(frame,dir)
    l = cos(el)*sin(az)
    m = cos(el)*cos(az)
    l,m
end
getlm(frame::ReferenceFrame,ra,dec) = getlm(frame,Direction("J2000",ra,dec))
getlm(frame::ReferenceFrame,source::Source) = getlm(frame,source.dir)

@doc """
Returns true if the source is above the horizon, false if the source
is below the horizon.
""" ->
function isabovehorizon(frame::ReferenceFrame,source::Source)
    az,el = getazel(frame,source)
    ifelse(el > 0.0Radian,true,false)
end

@doc """
Returns a list of sources that are above the horizon for the given
`ReferenceFrame`.
""" ->
function getsources(frame::ReferenceFrame)
    # TODO: Make this significantly better
    sources = Source[]
    push!(sources,Source("Cyg A",Direction("J2000",ra"19h59m17.24s",dec"+40d44m23.35s"),21850.24,47e6,[-0.51,-0.18]))
    push!(sources,Source("Cas A",Direction("J2000",ra"23h23m45.55s",dec"+58d25m25.39s"), 8628.39,47e6,[-0.14,-1.51]))

    sources_abovehorizon = Source[]
    for source in sources
        if isabovehorizon(frame,source)
            push!(sources_abovehorizon,source)
        end
    end
    sources_abovehorizon
end

