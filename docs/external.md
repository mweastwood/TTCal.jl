<!---
This is an auto-generated file and should not be edited directly.
-->

## API

### ConstantBeam

```
ConstantBeam <: BeamModel
```

In this beam model, the Jones matrix is assumed to be unity in every direction.

### DiagonalJonesMatrix

```
immutable DiagonalJonesMatrix
```

This type represents a Jones matrix that is diagonal.

These matrices are used to represent the complex gains of each antenna without accounting for the off-diagonal polarization leakage terms.

### GainCalibration

```
immutable GainCalibration <: Calibration
```

This type stores the information for calibrating the electronic gains of the interferometer. That is, it stores complex gains and flags for each antenna, frequency channel, and polarization.

```
GainCalibration(Nant, Nfreq)
```

Create a calibration table for `Nant` antennas with `Nfreq` frequency channels where all the gains are initially set to unity.

### HermitianJonesMatrix

```
immutable HermitianJonesMatrix
```

This type represents a Jones matrix that is Hermitian.

These matrices are useful for representing the $xx$, $xy$, $yx$, and $yy$ flux, because $xx$ and $yy$ are constrained to be real while $xy$ and $yx$ are complex conjugates.

### JonesMatrix

```
immutable JonesMatrix
```

This type represents a 2x2 complex Jones matrix.

\[ \begin{pmatrix}     j_{xx} & j_{xy} \\
    j_{yx} & j_{yy} \\
\end{pmatrix} \]

### MeasurementSet

```
immutable MeasurementSet
```

This type is a wrapper around `CasaCore.Tables.Table` that is intended to simplify most of the common interactions between TTCal and measurement sets.

```
MeasurementSet(name)
```

Open the measurement set at the given location. Assorted quantities that are commonly used by TTCal are automatically loaded and stored in fields.

### Memo178Beam

```
Memo178Beam <: BeamModel
```

This beam is based on the parametric fit to EM simulations presented in LWA memo 178 by Jayce Dowell.

[http://www.faculty.ece.vt.edu/swe/lwa/memo/lwa0178a.pdf]

### MuellerMatrix

```
immutable MuellerMatrix
```

This type represents a Mueller matrix.

```
MuellerMatrix(J::JonesMatrix)
```

Create a Mueller matrix from the given Jones matrix.

### Point



### PolarizationCalibration

```
immutable PolarizationCalibration <: Calibration
```

This type stores the information for calibrating the polarization of the interferometer. That is, it stores Jones matrices and flags for each antenna and each frequency channel.

```
PolarizationCalibration(Nant, Nfreq)
```

Create a calibration table for `Nant` antennas with `Nfreq` frequency channels where all the Jones matrices are initially set to the identity matrix.

### SineBeam

```
SineBeam <: BeamModel
```

This beam is azimuthally symmetric and independent of frequency. The gain of an individual dipole scales as $\sin(elevation)^\alpha$.

### Source

```
immutable Source
```

This type represents a radio source.

Each source is composed of one or more components.

### Spectrum

```
immutable Spectrum
```

These sources have a multi-component power-law spectrum such that:

\[     \log_{10} S = \log_{10} S_0 + \sum_{n=1}^N \alpha_n \log_{10}\left(\frac{\nu}{\nu_0}\right)^n \]

where $S$ is a Stokes parameter, and $\alpha_n$ is the list of spectral indices. At least one spectral index needs to be provided.

### StokesVector

```
immutable StokesVector
```

This type represents a Stokes vector.

\[ \begin{pmatrix}     I \\
    Q \\
    U \\
    V \\
\end{pmatrix} \]

### applycal!

```
applycal!(ms::MeasurementSet, calibration::Calibration;
          apply_to_corrected = false, force_imaging_columns = false)
```

Apply the calibration to the given measurement set.

**Arguments:**

  * `ms` - the measurement set to which the calibration will be applied
  * `calibration` - the calibration that will be applied

**Keyword Arguments:**

  * `apply_to_corrected` - if this is set to true, the calibration will be     applied to the CORRECTED_DATA column instead of the DATA column
  * `force_imaging_columns` - if this is set to true, the calibrated data     will be written to the CORRECTED_DATA column regardless of whether     or not the column already exists

### corrupt!

```
corrupt!(data::Array{Complex64,3}, cal::Calibration, ant1, ant2)
```

Corrupt the model data as if it had been observed with an instrument with the given calibration.

```
corrupt!(data::Array{Complex64,3}, flags::Array{Bool,3},
         cal::Calibration, ant1, ant2)
```

Corrupt the data as if it was observed with the given calibration.

### fitvis

```
fitvis(ms::MeasurementSet, direction::Direction;
       maxiter = 20, tolerance = 1e-3, minuvw = 0.0) -> l,m
```

Fit for the location of a point source near the given direction.

### gaincal

```
gaincal(ms::MeasurementSet, sources::Vector{Source}, beam::BeamModel;
        maxiter = 20, tolerance = 1e-3, minuvw = 0.0,
        reference_antenna = "1x", force_imaging_columns = false)
```

Solve for the interferometer's electronic gains.

**Arguments:**

  * `ms` - the measurement set from which to derive the calibration
  * `sources` - the list of points sources to use as the sky model
  * `beam` - the beam model

**Keyword Arguments:**

  * `maxiter` - the maximum number of Runge-Kutta steps to take on each     frequency channel
  * `tolerance` - the relative tolerance to use while checking to see if     more iterations are required
  * `minuvw` - the minimum baseline length (measured in wavelengths) to be     used during the calibration procedure
  * `reference_antenna` - a string containing the antenna number and polarization     whose phase will be chosen to be zero (eg. "14y" or "62x")
  * `force_imaging_columns` - if this is set to true, the MODEL_DATA column     will be created and populated with model visibilities even if it     doesn't already exist

### genvis

```
genvis(ms::MeasurementSet, sources, beam::BeamModel)
```

Generate model visibilities for the given list of sources and the given beam model.

No gridding is performed, so the runtime of this naive algorithm scales as $O(N_{base} \times N_{source})$.

### getspec

```
getspec(data, flags, l, m, u, v, w, ν, ant1, ant2) -> Vector{HermitianJonesMatrix}
```

Compute the spectrum of a source located at $(l,m)$ in all of the polarized correlation products.

```
getspec(ms::MeasurementSet, direction::Direction;
        minuvw = 0.0) -> Vector{HermitianJonesMatrix}
```

This function extracts the spectrum in a given direction by means of an inverse discrete Fourier transform.

Note that no gridding is performed, so this does *not* use a fast Fourier transform. However, the inverse discrete Fourier transform *is* the least squares estimator for the flux in a given direction (if all baselines are weighted equally).

### peel!

```
peel!{T<:Calibration}(::Type{T},
                      ms::MeasurementSet,
                      sources::Vector{Source},
                      beam::BeamModel;
                      peeliter = 3,
                      maxiter = 20,
                      tolerance = 1e-3,
                      minuvw = 0.0)
```

Peel the given list of sources from the measurement set.

The type supplied as the first argument determines the manner in which the sources are peeled:

  * `PolarizationCalibration` - each source receives a full set of Jones matrices
  * `GainCalibration` - each source receives a full set of complex gains
  * `AmplitudeCalibration` - each source receives a full set of gain amplitudes

### polcal

```
polcal(ms::MeasurementSet, sources::Vector{Source}, beam::BeamModel;
       maxiter = 20, tolerance = 1e-3, minuvw = 0.0,
       force_imaging_columns = false)
```

Solve for the polarization properties of the interferometer.

**Arguments:**

  * `ms` - the measurement set from which to derive the calibration
  * `sources` - the list of points sources to use as the sky model
  * `beam` - the beam model

**Keyword Arguments:**

  * `maxiter` - the maximum number of Runge-Kutta steps to take on each     frequency channel
  * `tolerance` - the relative tolerance to use while checking to see if     more iterations are required
  * `minuvw` - the minimum baseline length (measured in wavelengths) to be     used during the calibration procedure
  * `force_imaging_columns` - if this is set to true, the MODEL_DATA column     will be created and populated with model visibilities even if it     doesn't already exist

### readsources

```
readsources(filename) -> Vector{Source}
```

Read the list of point sources from the given JSON file. The format must be as follows:

```
[
    {
        "ref": "Baars et al. 1977",
        "name": "Cas A",
        "ra": "23h23m24s",
        "dec": "58d48m54s",
        "I": 555904.26,
        "freq": 1.0e6,
        "index": [-0.770]
    },
    {
        "ref": "Baars et al. 1977",
        "name": "Cyg A",
        "ra": "19h59m28.35663s",
        "dec": "+40d44m02.0970s",
        "I": 49545.02,
        "freq": 1.0e6,
        "index": [+0.085,-0.178]
    }
]
```

Additional Stokes parameters may also be specified.

```
{
    ...
    "I": 100
    "Q": 20
    "U": 3.14
    "V": -30
    ...
}
```

A right ascension and declination does not need to be specified if the name of the source is "Sun", "Moon", or "Jupiter". These sources will have their location automatically determined by CasaCore.

### subsrc!

```
subsrc!(ms::MeasurementSet, dir::Direction; minuvw = 0.0)
```

Subtract all of the measured flux from a given direction.

This can be used to remove RFI sources provided they have a known direction.

```
subsrc!(ms::MeasurementSet, sources::Vector{Source}, beam::BeamModel)
```

Remove the list of sources from the measurement set.

### writesources

```
writesources(filename, sources::Vector{Source})
```

Write the list of sources to the given location as a JSON file. These sources can be read back in again using the `readsources` function.

