module TTCal

export Interferometer, Source
export flagdata!
export bandpass
export fitvis
export applycal

using Devectorize
using JSON, NPZ
using SIUnits
using CasaCore.Measures
using CasaCore.Tables

include("units.jl")
include("rungekutta.jl")
include("interferometer.jl")
include("sourcemodel.jl")
include("visibilities.jl")

include("flagdata.jl")
include("bandpass.jl")
include("selfcal.jl")
include("applycal.jl")

#=
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
        sources  = getsources(frame)
        criteria = StoppingCriteria(args["niter"],
                                    args["tol"])
        gains = bandpass(lwa,ms,sources,criteria)
        save(args["gaintable"],gains)
    end
    nothing
end
=#

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

