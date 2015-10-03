<!---
This is an auto-generated file and should not be edited directly.
-->

## AmplitudeCalibration

```
AmplitudeCalibration <: Calibration
```

This type stores the information for calibrating only the amplitude of the electronic gains. It stores the gain amplitudes and flags for each antenna, frequency channel, and polarization.
## GainCalibration

```
GainCalibration <: Calibration
```

This type stores the information for calibrating the electronic gains of the interferometer. That is, it stores complex gains and flags for each antenna, frequency channel, and polarization.
## PointSource

These sources have a multi-component power-law spectrum such that:

```
log(flux) = log(I) + index[1]*log(ν/reffreq)
                   + index[2]*log²(ν/reffreq) + ...
```

Polarized fluxes are obtained in a similar manner by substituting Q/U/V for I in the above expression.
## PolarizationCalibration

```
PolarizationCalibration <: Calibration
```

This type stores the information for calibrating the polarization of the interferometer. That is, it stores Jones matrices and flags for each antenna and each frequency channel.
## ampcal

```
ampcal(ms::Table, sources::Vector{PointSource};
       maxiter = 30, tolerance = 1e-3,
       force_imaging_columns = false)
```

Solve for the amplitude of the interferometer's gains.
## applycal!

Apply the calibration to the given measurement set.
## corrupt!

Corrupt the model data as if it had been observed with an instrument with the given calibration.
## fitvis

Fit the visibilities to a model of point sources. The input model needs to have the positions of the point sources relatively close, but the flux can be wildly off.
## gaincal

```
gaincal(ms::Table, sources::Vector{PointSource};
        maxiter = 20, tolerance = 1e-5,
        force_imaging_columns = false,
        reference_antenna = 1)
```

Solve for the interferometer's electronic gains.
## genvis

Generate model visibilities for a given point source model. No gridding is performed, so the runtime of this naive algorithm scales as O(Nbase×Nsource).
## getspec

This function extracts the spectrum in a given direction by means of an inverse discrete Fourier transform. Note that no gridding is performed, so this does *not* use a fast Fourier transform. However, the inverse discrete Fourier transform *is* the least squares estimator for the flux in a given direction (if all baselines are weighted equally).
## peel!

```
peel!{T<:Calibration}(::Type{T},
                      ms::Table,
                      sources::Vector{PointSource};
                      maxiter = 20, tolerance = 1e-3,
                      minuvw = 0.0)
```

Peel the given list of sources from the measurement set.

The type supplied as the first argument determines the manner in which the sources are peeled:

  * `PolarizationCalibration` - each source receives a full set of Jones matrices
  * `GainCalibration` - each source receives a full set of complex gains
  * `AmplitudeCalibration` - each source receives a full set of gain amplitudes
## polcal

```
polcal(ms::Table, sources::Vector{PointSources};
       maxiter = 20, tolerance = 1e-5,
       force_imaging_columns = false)
```

Solve for the polarization properties of the interferometer.
## readsources


## subsrc!

```
subsrc!(ms::Table,dir::Direction)
```

Subtract all of the measured flux from a given direction.

This can be used to remove RFI sources provided they have a known direction.

```
subsrc!(ms::Table,sources::Vector{PointSource}
```

Remove the list of sources from the measurement set.
## writesources


