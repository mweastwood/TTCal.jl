# These are helper routines for interfacing with measurement sets

function uvw(ms::Table)
    uvw = ms["UVW"]
    u = addunits(uvw[1,:],Meter)
    v = addunits(uvw[2,:],Meter)
    w = addunits(uvw[3,:],Meter)
    u,v,w
end

function freq(ms::Table)
    spw = Table(ms[kw"SPECTRAL_WINDOW"])
    ν = addunits(spw["CHAN_FREQ",1],Hertz)
    ν
end

function reference_frame(ms::Table)
    frame = ReferenceFrame()
    set!(frame,Epoch("UTC",ms["TIME",1]*Second))
    set!(frame,Measures.observatory(frame,"OVRO_MMA"))
    frame
end

