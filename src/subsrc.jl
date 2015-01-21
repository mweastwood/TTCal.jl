################################################################################
# Public Interface

#using PyPlot
function subsrc!{T<:Union(Source,Vector{Source})}(ms::Table,sources::T)
    frame = ReferenceFrame()
    set!(frame,Epoch("UTC",ms["TIME",1]*Second))
    set!(frame,Measures.observatory(frame,"OVRO_MMA"))

    uvw  = ms["UVW"]
    u = addunits(squeeze(uvw[1,:],1),Meter)
    v = addunits(squeeze(uvw[2,:],1),Meter)
    w = addunits(squeeze(uvw[3,:],1),Meter)

    spw = Table(ms[kw"SPECTRAL_WINDOW"])
    ν = addunits(spw["CHAN_FREQ",1],Hertz)

    model = visibilities(frame,sources,u,v,w,ν)
    data  = ms["CORRECTED_DATA"]

    #figure(1); clf()
    #b = stripunits(sqrt(u.^2+v.^2+w.^2))
    #plot(b,squeeze(data[1,1,:],(1,2)),"ro")
    #plot(b,squeeze(model[1,1,:],(1,2)),"bo")

    ms["CORRECTED_DATA"] = data - model
end

