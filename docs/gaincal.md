## Gain Calibration

With this calibration TTCal will solve for one complex gain (amplitude and phase)
per antenna polarization per frequency channel. That is for complex $g_i$ and
$i$, $j$ labeling the pair of antennas composing the baseline, TTCal seeks to minimize

\[
    \sum_{i,j} \left\| V_{i,j} - g_i g_j^* M_{i,j} \right\|^2
\]

where $V_{i,j}$ is the measured visibility and $M_{i,j}$ is the model visibility
measured on the baseline $i$, $j$.

TTCal uses the iterative routine independently described by
[Mitchell et al. 2008](http://adsabs.harvard.edu/abs/2008ISTSP...2..707M) and
[Salvini & Wijnholds 2014](http://adsabs.harvard.edu/abs/2014A%26A...571A..97S)
to solve for the complex gains.

For more information on defining a sky model see [Sky Models](sources.md).

TTCal will flag frequency channels that do not converge within the specified maximum
number of iterations. However, TTCal will not attempt to identify bad antennas
that may be poisoning the calibration. These antennas must be flagged ahead of
time.

*TTCal assumes that the measurement set contains all baselines, multiple
frequency channels, but only a single integration.*

### Running from a Julia script

```julia
using TTCal
ms      = MeasurementSet("data.ms")
sources = readsources("sources.json")
beam    = TTCal.SineBeam()
cal = gaincal(ms, sources, beam, maxiter = 20, tolerance = 1e-3)
applycal!(ms, cal, force_imaging_columns = true)
unlock(ms)
```

### Running from the command line

For a list of all available options, run:
```
ttcal.jl gaincal --help
```

To calibrate a standard OVRO LWA dataset:
```
ttcal.jl gaincal --input data.ms --output calibration.jld --sources sources.json
ttcal.jl applycal --input data.ms --calibration calibration.jld
```

