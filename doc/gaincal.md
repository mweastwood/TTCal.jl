## Gain Calibration

With this calibration, TTCal will solve for one complex gain per antenna polarization
per frequency channel. TTCal assumes that the measurement set contains all baselines
for a single integration (multiple frequency channels are allowed though).
TTCal uses the iterative routine described by Mitchell et al. 2008 and
Salvini & Wijnholds 2014 to solve for the complex gains.

For example, to calibrate a standard OVRO LWA dataset:
```
ttcal.jl gaincal --input data.ms --output calibration.tt --sources sources.json
ttcal.jl applycal --input data.ms --calibration calibration.ttb
```
See the "Sky Models" section for an example `sources.json` file.
For a list of all available options, run:
```
ttcal.jl gaincal
```

TTCal will flag frequency channels that do converge within the specified maximum
number of iterations. However, TTCal will not attempt to identify bad antennas
that may be poisoning the calibration. These antennas must be flagged ahead of
time.

