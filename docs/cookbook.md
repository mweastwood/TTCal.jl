# Cookbook

## Gain Calibration

**Running from a Julia script**
``` julia
using CasaCore.Tables
using TTCal

# open the measurement set and read the DATA column
ms = Table("data.ms")
visibilities = TTCal.read(ms, "DATA")

# read metadata from the measurement set and choose a beam model
meta = Metadata(ms)
beam = SineBeam()

# read the sky model
sources = readsources("sources.json")

# solve for the calibration
cal = gaincal(visibilities, meta, beam, sources)

# apply the calibration to the visibilities
applycal!(visibilities, meta, cal)

# write the corrected visibilities to the CORRECTED_DATA column
TTCal.write(ms, "CORRECTED_DATA", visibilities)

# release the lock on the measurement set so that other processes can use it
unlock(ms)
```

**Running from the command line**
``` sh
ttcal.jl gaincal --input data.ms --output calibration.jld --sources sources.json
ttcal.jl applycal --input data.ms --calibration calibration.jld
```

