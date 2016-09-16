# Cookbook

## Introduction

Most of these Julia language examples below will require some amount of setup. Use the following
code to do the necessary setup.

``` julia
using CasaCore.Tables
using TTCal

# open the measurement set and read the DATA column
ms = Table("data.ms")
visibilities = TTCal.read(ms, "DATA")

# read metadata from the measurement set and choose a beam model
meta = Metadata(ms) # contains information about antenna positions, frequency channels, etc.
beam = SineBeam()

# read the sky model
sources = readsources("sources.json")
```

Additionally when opening a measurement set we acquire a lock that prevents other processes from
interacting with it until we have finished. This lock will automatically be released when the julia
process exits, but it is generally good practice to release the lock yourself at the end of the
script.

``` julia
unlock(ms)
```

## Gain Calibration

**Running from a Julia script**
``` julia
# ... setup ...
# solve for the calibration
cal = gaincal(visibilities, meta, beam, sources)
# apply the calibration to the visibilities
applycal!(visibilities, meta, cal)
# write the corrected visibilities to the CORRECTED_DATA column
TTCal.write(ms, "CORRECTED_DATA", visibilities)
```

**Running from the command line**
``` sh
ttcal.jl gaincal data.ms calibration.jld sources.json
ttcal.jl applycal data.ms calibration.jld
```

## Polarization Calibration

**Running from a Julia script**
``` julia
# ... setup ...
# solve for the calibration
cal = polcal(visibilities, meta, beam, sources)
# apply the calibration to the visibilities
applycal!(visibilities, meta, cal)
# write the corrected visibilities to the CORRECTED_DATA column
TTCal.write(ms, "CORRECTED_DATA", visibilities)
```

**Running from the command line**
``` sh
ttcal.jl polcal data.ms calibration.jld sources.json
ttcal.jl applycal data.ms calibration.jld
```

## Peeling

**Running from the command line**
``` sh
ttcal.jl peel data.ms sources.json
```

