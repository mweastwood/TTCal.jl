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

## Command Line Interface

    $ ttcal.jl flagdata --input test.ms --antennas 8 16 24

This will flag all the baselines containing antennas 8, 16, and 24 (where antenna 1 is the first antenna).

    $ ttcal.jl bandpass --input test.ms --output gains.npy --sources sources.json

This will derive gains for the given measurement set using a model sky specified by the point sources in sources.json. The gains are stored as numpy arrays, which should make it easy to read and manipulate these gains from within python. An example JSON file specifying the sky model is included at the bottom of this email.

    $ ttcal.jl applycal --input test.ms --calibration gains.npy

This will apply the calibration derived from the previous example. At the moment, if a CORRECTED_DATA column exists, it will be used, otherwise the DATA column will be overwritten with the corrected visibilities.


```
[
    {
        "ref": "Baars et al. 1977",
        "name": "Cas A",
        "ra": "23h23m24s",
        "dec": "58d48m54s",
        "flux": 555904.26,
        "freq": 1.0e6,
        "index": [-0.770]
    }, 
    {
        "ref": "Baars et al. 1977",
        "name": "Cyg A",
        "ra": "19h59m28.35663s",
        "dec": "+40d44m02.0970s",
        "flux": 49545.02,
        "freq": 1.0e6,
        "index": [+0.085,-0.178]
    }
]
```
