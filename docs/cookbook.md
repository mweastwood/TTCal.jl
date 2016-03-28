# Cookbook

## Gain Calibration

**Running from a Julia script**
``` julia
using CasaCore.Tables
using TTCal

# open the measurement set and read the DATA column
ms = Table("data.ms")
visibilities = TTCal.get_data(ms)

# read metadata from the measurement set and choose a beam model
beam = SineBeam()
meta = collect_metadata(ms, beam)

# read the sky model
sources = readsources("sources.json")

# solve for the calibration
cal = gaincal(visibilities, meta, sources)

# apply the calibration to the visibilities
applycal!(visibilities, meta, cal)

# write the corrected visibilities to the CORRECTED_DATA column
TTCal.set_corrected_data!(ms, visibilities)

# release the lock on the measurement set so that other processes can use it
unlock(ms)
```

**Running from the command line**
``` sh
ttcal.jl gaincal --input data.ms --output calibration.jld --sources sources.json
ttcal.jl applycal --input data.ms --calibration calibration.jld
```

