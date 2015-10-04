## Polarization Calibration

For a polarization calibration, TTCal solves for one Jones matrix per antenna per
frequency channel. This is accomplished by using the fully polarized version of
the iterative routine described by Mitchell et al. 2008 and Salvini & Wijnholds 2014.

A standard OVRO LWA dataset can be calibrated with:
```
ttcal.jl polcal --input data.ms --output calibration.tt --sources sources.json
ttcal.jl applycal --input data.ms --calibration calibration.tt
```

