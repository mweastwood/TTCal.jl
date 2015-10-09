# TTCal

![TTCal](ttcal.png)

TTCal is a calibration routine developed for the OVRO LWA.

## Getting Started

TTCal requires the latest development version of the [Julia](http://julialang.org/) programming language.
```
$ git clone https://github.com/JuliaLang/julia.git
$ cd julia
$ make -j
```

To install TTCal, run:
```
$ julia
julia> Pkg.clone("CasaCore")
julia> Pkg.build("CasaCore")
julia> Pkg.clone("https://github.com/mweastwood/TTCal.jl.git")
julia> Pkg.test("TTCal")
```
If all the tests pass, you are ready to begin using TTCal.
Simply add the `ttcal.jl` file to your `PATH` environment variable.
You can see the list of available commands by running:
```
$ ttcal.jl --help
usage: ttcal.jl [-h] {gaincal|polcal|peel|applycal}

A calibration routine developed for the OVRO LWA.

commands:
  gaincal     Solve for a gain calibration.
  polcal      Solve for a polarization calibration.
  peel        Peel sources from the dataset.
  applycal    Apply a calibration.

optional arguments:
  -h, --help  show this help message and exit
```

