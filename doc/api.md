<!---
This is an auto-generated file and should not be edited directly.
-->

## AmplitudeCalibration

```
AmplitudeCalibration <: Calibration
```

This type stores the information for calibrating only the amplitude of the electronic gains. It stores the gain amplitudes and flags for each antenna, frequency channel, and polarization.

```
AmplitudeCalibration(Nant,Nfreq)
```

Create a calibration table for `Nant` antennas with `Nfreq` frequency channels where all the amplitudes are initially set to unity.

## GainCalibration

```
GainCalibration <: Calibration
```

This type stores the information for calibrating the electronic gains of the interferometer. That is, it stores complex gains and flags for each antenna, frequency channel, and polarization.

```
GainCalibration(Nant,Nfreq)
```

Create a calibration table for `Nant` antennas with `Nfreq` frequency channels where all the gains are initially set to unity.

## MeasurementSet

```
immutable MeasurementSet
```

This type is a wrapper around `CasaCore.Tables.Table` that is intended to simplify most of the common interactions between TTCal and measurement sets.

```
MeasurementSet(name)
```

Open the measurement set at the given location. Assorted quantities that are commonly used by TTCal are automatically loaded and stored in fields.

## PointSource

```
type PointSource
```

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

```
PolarizationCalibration(Nant,Nfreq)
```

Create a calibration table for `Nant` antennas with `Nfreq` frequency channels where all the Jones matrices are initially set to the identity matrix.

## ampcal

```
ampcal(ms::MeasurementSet,
       sources::Vector{PointSource},
       beam::BeamModel;
       maxiter = 30, tolerance = 1e-3,
       force_imaging_columns = false)
```

Solve for the amplitude of the interferometer's gains.

## applycal!

Apply the calibration to the given measurement set.

## corrupt!

Corrupt the model data as if it had been observed with an instrument with the given calibration.

## fitvis

```
fitvis(ms::MeasurementSet,
       sources::Vector{PointSource};
       maxiter::Int = 20,
       tolerance::Float64 = 1e-3,
       minuvw::Float64 = 0.0)
```

Fit for the location of each point source.

## gaincal

```
gaincal(ms::MeasurementSet,
        sources::Vector{PointSource},
        beam::BeamModel;
        maxiter = 20, tolerance = 1e-5,
        force_imaging_columns = false,
        reference_antenna = 1)
```

Solve for the interferometer's electronic gains.

## genvis

```
genvis(ms::MeasurementSet,
       sources::Union{PointSource,Vector{PointSource}},
       beam::BeamModel)
```

Generate model visibilities for the given list of sources and the given beam model.

No gridding is performed, so the runtime of this naive algorithm scales as $O(N_{base} \times N_{source})$.

## getspec

```
getspec(ms::MeasurementSet, dir::Direction) -> xx,xy,yx,yy
```

This function extracts the spectrum in a given direction by means of an inverse discrete Fourier transform.

Note that no gridding is performed, so this does *not* use a fast Fourier transform. However, the inverse discrete Fourier transform *is* the least squares estimator for the flux in a given direction (if all baselines are weighted equally).

```
getspec(data, flags, l, m, u, v, w, ν, ant1, ant2) -> xx,xy,yx,yy
```

Compute the spectrum of a source located at $(l,m)$ in all of the polarized correlation products.

## peel!

```
peel!{T<:Calibration}(::Type{T},
                      ms::MeasurementSet,
                      sources::Vector{PointSource},
                      beam::BeamModel;
                      maxiter = 20,
                      tolerance = 1e-3,
                      minuvw = 0.0)
```

Peel the given list of sources from the measurement set.

The type supplied as the first argument determines the manner in which the sources are peeled:

  * `PolarizationCalibration` - each source receives a full set of Jones matrices
  * `GainCalibration` - each source receives a full set of complex gains
  * `AmplitudeCalibration` - each source receives a full set of gain amplitudes

## polcal

```
polcal(ms::MeasurementSet,
       sources::Vector{PointSources},
       beam::BeamModel;
       maxiter = 20, tolerance = 1e-5,
       force_imaging_columns = false)
```

Solve for the polarization properties of the interferometer.

## readsources



## subsrc!

```
subsrc!(ms::MeasurementSet,
        sources::Vector{PointSource},
        beam::BeamModel)
```

Remove the list of sources from the measurement set.

```
subsrc!(ms::MeasurementSet, dir::Direction)
```

Subtract all of the measured flux from a given direction.

This can be used to remove RFI sources provided they have a known direction.

## writesources



