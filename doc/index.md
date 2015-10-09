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
julia> Pkg.clone("https://github.com/mweastwood/CLI.jl.git")
julia> Pkg.clone("https://github.com/mweastwood/CasaCore.jl.git")
julia> Pkg.clone("https://github.com/mweastwood/TTCal.jl.git")
julia> Pkg.test("TTCal")
```
If all the tests pass, you are ready to begin using TTCal.
Simply add the `ttcal.jl` file to your `PATH` environment variable.
You can see the list of available commands by running:
```
$ ttcal.jl 
 TTCal
=======
A calibration routine developed for the OVRO LWA.
Written by Michael Eastwood (mweastwood@astro.caltech.edu).

usage: ttcal.jl command options...

commands:
  gaincal         Solve for a gain calibration.
  polcal          Solve for a polarization calibration.
  peel            Peel sources from the dataset.
  applycal        Apply a calibration.

Please provide one of the listed commands.
```

