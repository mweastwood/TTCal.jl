module TTCal

import Base: show
using CasaCore.Measures
using CasaCore.Tables
using SIUnits

include("units.jl")
include("sourcemodel.jl")
include("interferometer.jl")
include("visibilities.jl")
include("bandpass.jl")
include("applycal.jl")
include("selfcal.jl")

function run(args)
    ms = [Table(ASCIIString(input)) for input in args["measurementsets"]]

    Nbase = numrows(ms[1])
    Nant  = round(Int,(sqrt(8Nbase+1)-1)/2)
    spw = Table(ms[1][kw"SPECTRAL_WINDOW"])
    Nfreq = length(spw["CHAN_FREQ",1])
    lwa = Interferometer(Nant,Nfreq,args["refant"],args["flags"])

    applycal_flag = args["applycal"]
    if applycal_flag
        gains = load(args["gaintable"])
        applycal(lwa,ms,gains)
    else
        frame = ReferenceFrame()
        set!(frame,Epoch("UTC",ms[1]["TIME",1]*Second))
        set!(frame,Measures.observatory(frame,"OVRO_MMA"))
        sources = getsources(frame)
        options = BandpassOptions(args["niter"],
                                  args["tol"],
                                  args["RK"],
                                  args["doubleprecision"])
        gains = bandpass(lwa,ms,sources,options)
        save(args["gaintable"],gains)
    end
    nothing
end

using NPZ
@doc """
Save the gains as a numpy array.
""" ->
function save(filename,gains)
    npzwrite(filename,gains)
    nothing
end
@doc """
Load the gains from a numpy array.
""" ->
function load(filename)
    npzread(filename)
end

end

