# These are helper routines for interfacing with measurement sets

function uvw(ms::Table)
    uvw_arr = ms["UVW"]
    u = addunits(squeeze(uvw_arr[1,:],1),Meter)
    v = addunits(squeeze(uvw_arr[2,:],1),Meter)
    w = addunits(squeeze(uvw_arr[3,:],1),Meter)
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
    set!(frame,observatory(frame,"OVRO_MMA"))
    frame
end

function ants(ms::Table)
    ant1 = ms["ANTENNA1"]
    ant2 = ms["ANTENNA2"]
    # (the +1 converts to a 1-based indexing scheme)
    ant1+1,ant2+1
end

