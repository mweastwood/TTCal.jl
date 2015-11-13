<!---
This is an auto-generated file and should not be edited directly.
-->

## Internal Documentation

### RK{N}

This singleton type is used to indicate which Runge-Kutta method should be used. For example, `RK{4}` tells us to use the RK4 method.

### congruence_transform

```
congruence_transform(J::JonesMatrix, K::HermitianJonesMatrix) -> J*K*J'
```

Compute the congruence transformation of $K$ with respect to $J$:

\[     K \rightarrow JKJ^* \]

### direction_cosines

```
direction_cosines(phase_dir::Direction{dir"J2000"}, dir::Direction{dir"J2000"}) -> l,m
```

Compute the direction cosines $(l,m)$ for the given direction with respect to the phase direction.

### fixphase!

```
fixphase!(cal::Calibration, reference_antenna)
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

### fringepattern

```
fringepattern(ϕ, Δϕ, N)
```

Compute $\exp(i(\phi+n\Delta\phi))$ where $\phi$, $\Delta\phi$, and $n = 1,\ldots,N$ define an equally spaced grid of points.

Using the sine and cosine angle addition rules, you can define an iterative method such that you only need to compute sines and cosines for a single iteration.

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
invert(cal::Calibration)
```

Returns the inverse of the given calibration. The Jones matrix $J$ of each antenna is set to $J^{-1}$.

### j2000_radec

```
j2000_radec(frame::ReferenceFrame, component::Component) -> ra,dec
```

Compute the J2000 right ascension and declination of the component (in radians).

### latlong_to_utm

```
latlong_to_utm(zone, latitude, longitude)
```

Convert the latitude and longitude (both in degrees) to an easting and a northing (both in meters).

This function is entirely based on the [Universal Transverse Mercator coordinate system](https://en.wikipedia.org/w/index.php?title=Universal_Transverse_Mercator_coordinate_system&oldid=683193579#Simplified_formulas) Wikipedia page.

### linear

```
linear(stokes::StokesVector) -> HermitianJonesMatrix
```

Take the set of Stokes parameters $(I,Q,U,V)$ and convert it to the set of correlations $(xx,xy,yx,yy)$.

### local_azel

```
local_azel(frame::ReferenceFrame, component::Component) -> az,el
```

Compute the local azimuth and elevation of the component (in radians).

### makesquare

```
makesquare(data, flags, ant1, ant2)
```

Pack the data into a square Hermitian matrix such that the data is ordered as follows:

\[     \begin{pmatrix}         V_{11} & V_{12} & V_{13} &        & \\
        V_{21} & V_{22} & V_{23} &        & \\
        V_{31} & V_{32} & V_{33} &        & \\
               &        &        & \ddots & \\
    \end{pmatrix} \]

Flagged correlations and autocorrelations are set to zero.

### rkstep

Take a Runge-Kutta step. * `step(x,args...)` must return the list of steps to take from the given location `x` * `rk` is the order of the Runge-Kutta method to use (eg. `RK4`) * `x` is the starting location * `args` is simply passed on as the second argument to `step`

### stefcal_step

```
stefcal_step(input,data,model) -> step
```

Given the `data` and `model` visibilities, and the current guess for the Jones matrices, solve for `step` such that the new value of the Jones matrices is `input+step`.

The update step is defined such that the new value of the Jones matrices minimizes

\[     \sum_{i,j}\|V_{i,j} - G_i M_{i,j} G_{j,new}^*\|^2$, \]

where $i$ and $j$ label the antennas, $V$ labels the measured visibilities, $M$ labels the model visibilities, and $G$ labels the Jones matrices.

*References:*

  * Michell, D. et al. 2008, JSTSP, 2, 5.
  * Salvini, S. & Wijnholds, S. 2014, A&A, 571, 97.

### stokes

```
stokes(correlations::HermitianJonesMatrix) -> StokesVector
```

Take the set of correlations $(xx,xy,yx,yy)$ and convert it to the set of Stokes parameters $(I,Q,U,V)$.

### utm_to_latlong

```
utm_to_latlong(zone, easting, northing,
               hemisphere = north)
```

Convert the easting and northing (both in meters) to a latitude and longitude (both in degrees).

This function is entirely based on the [Universal Transverse Mercator coordinate system](https://en.wikipedia.org/w/index.php?title=Universal_Transverse_Mercator_coordinate_system&oldid=683193579#Simplified_formulas) Wikipedia page.

