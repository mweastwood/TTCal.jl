<!---
This is an auto-generated file and should not be edited directly.
-->

## Internal Documentation

### ConstantBeam

```
ConstantBeam <: BeamModel
```

In this beam model, the Jones matrix is assumed to be unity in every direction.

### JonesMatrix

This type represents a 2x2 complex Jones matrix.

```
⌈xx xy⌉
⌊yx yy⌋
```

### Memo178Beam

```
Memo178Beam <: BeamModel
```

This beam is based on the parametric fit to EM simulations presented in LWA memo 178 by Jayce Dowell.

[http://www.faculty.ece.vt.edu/swe/lwa/memo/lwa0178a.pdf]

### RK{N}

This singleton type is used to indicate which Runge-Kutta method should be used. For example, `RK{4}` tells us to use the RK4 method.

### SineBeam

```
SineBeam <: BeamModel
```

This beam is azimuthally symmetric and independent of frequency. The gain of an individual dipole scales as $\sin(elevation)^\alpha$.

### ampcal_step

```
ampcal_step(amplitudes,data,model) -> step
```

Given the `data` and `model` visibilities, and the current guess for the gain `amplitudes`, solve for `step` such that the new value of the amplitudes is `amplitudes+step`.

### dir2azel

```
dir2azel(frame::ReferenceFrame, dir::Direction)
```

Convert the direction into a local azimuth and elevation.

### dir2lm

```
dir2lm(phase_dir::Direction{dir"J2000"}, dir::Direction{dir"J2000"})
```

Convert the direction into the standard radio coordinate system.

### dir2radec

```
dir2radec(frame::ReferenceFrame, dir::Direction)
```

Convert the direction into a J2000 right ascension and declination.

### fixphase!

```
fixphase!(cal::AmplitudeCalibration, reference_antenna)
```

With an amplitude calibration there is no freedom to pick an arbitrary phase, so this function does nothing.

```
fixphase!(cal::GainCalibration, reference_antenna)
```

Set the phase of the reference antenna and polarization to zero.

**Arguments:**

  * `cal` - the calibration that will have its phase adjusted
  * `reference_antenna` - a string containing the antenna number and polarization     whose phase will be chosen to be zero (eg. "14y" or "62x")

### flag_short_baselines!

```
flag_short_baselines!(flags, minuvw, u, v, w, ν)
```

Flag all of the baselines whose length is less than `minuvw` wavelengths.

This is a common operation that can mitigate contamination by unmodeled diffuse emission.

### force_to_horizon

```
force_to_horizon(l,m)
```

This function forces the coordinates $(l,m)$ to be above the horizon.

Although this is a nasty hack, it is necessary for fitting some sources that are near the horizon.

### fringepattern!

Compute exp(i(ϕ+nΔϕ)) where ϕ and Δϕ define an equally space grid of points where n = 1 to N.

Using the sine and cosine angle addition rules, you can define an iterative method such that you only need to compute sines and cosines for a single iteration.

### gaincal_step

```
gaincal_step(gains,data,model) -> step
```

Given the `data` and `model` visibilities, and the current guess for the electronic `gains`, solve for `step` such that the new value of the gains is `gains+step`.

The update step is defined such that the new value of the gains minimizes

[     \sum_{i,j}\|V_{i,j} - g_i g_{j,new}^* M_{i,j}\|^2$, ]

where $i$ and $j$ label the antennas, $V$ labels the measured visibilities, $M$ labels the model visibilities, and $g$ labels the complex gains.

*References:*

  * Michell, D. et al. 2008, JSTSP, 2, 5.
  * Salvini, S. & Wijnholds, S. 2014, A&A, 571, 97.

### get_corrected_data

```
get_corrected_data(ms::MeasurementSet)
```

Get the CORRECTED_DATA column if it exists. Otherwise settle for the DATA column.

### get_flags

```
get_flags(ms::MeasurementSet)
```

Get the flags from the dataset, but this information is stored in multiple locations. Unify all these flags before returning.

### invert

```
invert(cal::ScalarCalibration)
```

Returns the inverse of the given calibration. The gain $g$ of each antenna is set to $1/g$.

```
invert(cal::PolarizationCalibration)
```

Returns the inverse of the given calibration. The Jones matrix $J$ of each antenna is set to $J^{-1}$.

### latlong_to_utm

```
latlong_to_utm(zone, latitude, longitude)
```

Convert the latitude and longitude (both in degrees) to an easting and a northing (both in meters).

This function is entirely based on the [Universal Transverse Mercator coordinate system](https://en.wikipedia.org/w/index.php?title=Universal_Transverse_Mercator_coordinate_system&oldid=683193579#Simplified_formulas) Wikipedia page.

### linear

```
linear(stokes) -> [xx,xy,yx,yy]
```

Take the vector of Stokes parameters `[I,Q,U,V]` and convert it to the vector of linear correlations.

### mueller

```
mueller(J::JonesMatrix)
```

Create a Mueller matrix from the given Jones matrix.

```
mueller(J1::JonesMatrix, J2::JonesMatrix)
```

Create a Mueller matrix from the two Jones matrices.

### polcal_makesquare

```
polcal_makesquare(data, flags, ant1, ant2)
```

Pack the data into a square Hermitian matrix where each element is a Jones matrix.

Compare this to the regular complex gain calibration where each element is a complex scalar.

The packing order is:

```
11 12 13
21 22 23
31 32 33
         .
           .
             .
```

### polcal_step

```
polcal_step(jones,data,model) -> step
```

Given the `data` and `model` visibilities, and the current guess for the Jones matrices (`jones`), solve for `step` such that the new value of the Jones matrices is `jones+step`.

### rkstep

Take a Runge-Kutta step. * `step(x,args...)` must return the list of steps to take from the given location `x` * `rk` is the order of the Runge-Kutta method to use (eg. `RK4`) * `x` is the starting location * `args` is simply passed on as the second argument to `step`

### scalar_makesquare

```
scalar_makesquare(data, flags, ant1, ant2)
```

Pack the data into a square Hermitian matrix such that the data is ordered as follows:

```
x₁x₁ x₁y₁ x₁x₂ x₁y₂
y₁x₁ y₁y₁ y₁x₂ y₁y₂
x₂x₁ y₂y₁ y₂x₂ y₂y₂
y₂x₁ y₂y₁ y₂x₂ y₂y₂
                     .
                       .
                         .
```

Flagged correlations and autocorrelations are set to zero.

### stokes

```
stokes(correlations) -> [I,Q,U,V]
```

Take the vector of correlations `[xx,xy,yx,yy]` and convert it to the Stokes parameters.

### utm_to_latlong

```
utm_to_latlong(zone, easting, northing,
               hemisphere = north)
```

Convert the easting and northing (both in meters) to a latitude and longitude (both in degrees).

This function is entirely based on the [Universal Transverse Mercator coordinate system](https://en.wikipedia.org/w/index.php?title=Universal_Transverse_Mercator_coordinate_system&oldid=683193579#Simplified_formulas) Wikipedia page.

### write_ds9_regions

```
write_ds9_regions(filename, sources::Vector{PointSource})
```

Write the list of sources to a DS9 region file. This file can then be loaded into DS9 to highlight sources within images.

