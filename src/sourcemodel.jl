# The source flux is modeled as follows:
#     log(S) = log(flux) + index[1]*log(ν/reffreq)
#                        + index[2]*log²(ν/reffreq) + ...
type Source
    name::ASCIIString
    ra::Quantity{Float64,ASCIIString}
    dec::Quantity{Float64,ASCIIString}
    flux::Float64
    reffreq::Float64
    index::Vector{Float64}
end

Source(name,ra,dec,flux,reffreq,index::Float64) = Source(name,ra,dec,flux,reffreq,[index])

function getflux(source::Source,frequency::Float64)
    logflux = log10(source.flux)
    for (i,index) in enumerate(source.index)
        logflux += index*log10(frequency/source.reffreq).^i
    end
    10.0.^logflux
end

function getflux(source::Source,frequencies::Vector{Float64})
    [getflux(source,frequency) for frequency in frequencies]
end

@doc """
Returns true if the source is above the horizon, false if the source
is below the horizon. The location of the observatory is hardcoded
to the OVRO MMA, but the local time is read from the measurement set.
""" ->
function isabovehorizon(ms::MeasurementSet,source::Source)
    frame = ReferenceFrame()
    set!(frame,getTime(ms))
    set!(frame,observatory(frame,"OVRO_MMA"))
    dir  = Direction("J2000",source.ra,source.dec)
    azel = measure(frame,"AZEL",dir)
    el   = azel.m[2]
    ifelse(el.value > 0.,true,false)
end

@doc """
Returns a list of sources that are above the horizon for the given
measurement set.
""" ->
function getsources(ms::MeasurementSet)
    # TODO: Make this significantly better
    sources = Source[]
    push!(sources,Source("Cyg A",q"19h59m17.24s",q"+40d44m23.35s",21850.24,47e6,[-0.51,-0.18]))
    push!(sources,Source("Cas A",q"23h23m45.55s",q"+58d25m25.39s", 8628.39,47e6,[-0.14,-1.51]))
    push!(sources,Source("Her A",q"16h48m56.65s",q"+05d27m26.30s",  797.79,47e6,[-1.02,+0.08]))
    push!(sources,Source("Tau A",q"05h34m31.94s",q"+22d00m52.2s",1770.,80e6,-0.27))
    push!(sources,Source("Vir A",q"12h30m49.42s",q"+12d23m28.0s",2400.,80e6,-0.86))
    push!(sources,Source("SNR G078.2+02.1",q"20h20m23.73s",q"+40d16m41.34s",  613.85,47e6,[+1.90,-7.97]))
    push!(sources,Source("SNR G082.2+05.3",q"20h19m35.82s",q"+45d33m07.52s",  328.09,47e6,[+0.36,-7.07]))
    push!(sources,Source("SNR G089.0+04.7",q"20h46m00.16s",q"+50d37m13.69s",  434.71,47e6,[-1.21,-0.93]))

    sources_abovehorizon = Source[]
    for source in sources
        if isabovehorizon(ms,source)
            push!(sources_abovehorizon,source)
        end
    end
    sources_abovehorizon
end

function show(io::IO,source::Source)
    print(io,@sprintf("%20s    %15s    %15s    %10.2fJy    %5.2fMHz",
                      source.name,source.ra,source.dec,source.flux,source.reffreq/1e6))
    for index in source.index
        print(io,@sprintf("    %+5.2f",index))
    end
    nothing
end

function show(io::IO,sources::Vector{Source})
    for (i,source) in enumerate(sources)
        show(io,source)
        i < length(sources) && print(io,'\n')
    end
    nothing
end

