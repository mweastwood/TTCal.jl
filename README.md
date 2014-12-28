# TTCal

[![Build Status](https://travis-ci.org/mweastwood/TTCal.jl.svg?branch=master)](https://travis-ci.org/mweastwood/TTCal.jl)

TTCal is a bandpass calibration routine developed for the OVRO LWA. Its primary advantage is that it is faster and simpler than Casa's equivalent bandpass task.

## Getting Started

This package requires the (currently unregistered) `CasaCore` package. To install `CasaCore`, run:
```julia
Pkg.clone("https://github.com/mweastwood/CasaCore.jl.git")
Pkg.build("CasaCore")
```
To install `TTCal`, run:
```julia
Pkg.clone("https://github.com/mweastwood/CasaCore.jl.git")
Pkg.test("TTCal")
```
If all the tests pass, you are ready to begin using `TTCal`.

## Unpolarized Bandpass Calibration
### Algorithm

## Polarized Bandpass Calibration

Not yet implemented.

## Direction-Dependent Calibration

Not yet implemented.

## Self-Calibration
