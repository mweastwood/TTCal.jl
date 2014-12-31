module TTCal

import Base: show
using CasaCore.Measures
using CasaCore.Tables

include("sourcemodel.jl")
include("interferometer.jl")
include("visibilities.jl")
include("bandpass.jl")
include("applycal.jl")

function run(args)
    ms = [Table(ASCIIString(input)) for input in args["measurementsets"]]

    lwa = LWA()
    lwa.flaggedantennas = args["flags"]
    lwa.refant = args["refant"]

    applycal_flag = args["applycal"]
    if applycal_flag
        gains = load(args["gaintable"])
        applycal(lwa,ms,gains)
    else
        sources = getsources(ms[1])
        options = BandpassOptions(args["niter"],
                                  args["tol"],
                                  args["RK"],
                                  args["doubleprecision"])
        bandpass(lwa,ms,sources,options)
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

