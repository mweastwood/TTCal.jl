# TTCal

[![Build Status](https://travis-ci.org/mweastwood/TTCal.jl.svg?branch=master)](https://travis-ci.org/mweastwood/TTCal.jl)
[![Coverage Status](https://coveralls.io/repos/mweastwood/TTCal.jl/badge.svg?branch=master)](https://coveralls.io/r/mweastwood/TTCal.jl?branch=master)

TTCal is a bandpass calibration routine developed for the OVRO LWA. Its primary advantage is that it is faster and simpler than Casa's equivalent bandpass task.

## Getting Started

To install `TTCal`, run:
```julia
Pkg.clone("https://github.com/mweastwood/TTCal.jl.git")
Pkg.test("TTCal")
```
If all the tests pass, you are ready to begin using `TTCal`.

## Unpolarized Bandpass Calibration
### Algorithm
Given a set of measured visibilities `V_ij` and model visibilities `M_ij`, the relationship between `V_ij` and `M_ij` is simply:
```
V_ij = g_i conj(g_j) M_ij
```
where `g_i` is the complex gain associated with antenna `i`. Using the identity `conj(V_ij) = V_ji`, TTCal packs the measured and model visibilities into square Hermitian matrices. In the absence of noise, the principle eigenvector of `V_ij/M_ij` is the complex gain for each antenna. TTCal uses this principle eigenvector as the starting point for an iterative solver.

TTCal's iterative solver is inspired by Stefcal.

## Polarized Bandpass Calibration
### Algorithm

## Direction-Dependent Calibration

Not yet implemented.

## Self-Calibration
