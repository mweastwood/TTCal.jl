module TTCal

import Base: show
using CasaCore

include("sourcemodel.jl")
include("interferometer.jl")
include("visibilities.jl")
include("bandpass.jl")
include("applycal.jl")

function run(args)
    ms = [MeasurementSet(ASCIIString(input)) for input in args["measurementsets"]]
    println(typeof(ms))
    if args["applycal"]
        # Apply the calibration only
        gains = load(args["gaintable"])
        applycal(ms,gains)
    else
        # Array parameters
        flaggedantennas = readdlm("flags.dat")
        lwa = LWA()
        lwa.flaggedantennas = args["flags"]
        lwa.refant = args["refant"]

        # Calibrate
        sources = getsources(ms[1])
        println(sources)
        options = BandpassOptions(args["niter"],
                                  args["tol"],
                                  args["RK"],
                                  args["doubleprecision"])
        gains = bandpass(lwa,ms,sources,options)
        applycal(lwa,ms,gains)
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

