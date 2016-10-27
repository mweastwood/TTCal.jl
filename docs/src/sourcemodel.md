# Sky Models

## Point Sources

For the OVRO LWA the two most important point sources for calibration are
[Cas A](https://en.wikipedia.org/wiki/Cassiopeia_A) and
[Cyg A](https://en.wikipedia.org/wiki/Cygnus_A), which are the two brightest
sources on the sky. The optimal time to solve for a calibration is when both
of these sources are above the horizon.

Currently TTCal expects the sky model to be specified as a list of point sources.
This is accomplished by listing the sources in a [JSON](http://www.json.org/) file.
For example the following file defines a sky model consisting of Cas A and Cyg A
using the published spectra for [Baars et al. 1977](http://adsabs.harvard.edu/abs/1977A%26A....61...99B).

```
[
    {
        "ref": "Baars et al. 1977",
        "name": "Cas A",
        "ra": "23h23m24s",
        "dec": "58d48m54s",
        "I": 555904.26,
        "Q": 0.0,
        "U": 0.0,
        "V": 0.0,
        "freq": 1.0e6,
        "index": [-0.770]
    },
    {
        "ref": "Baars et al. 1977",
        "name": "Cyg A",
        "ra": "19h59m28.35663s",
        "dec": "+40d44m02.0970s",
        "I": 49545.02,
        "Q": 0.0,
        "U": 0.0,
        "V": 0.0,
        "freq": 1.0e6,
        "index": [+0.085,-0.178]
    }
]
```

There are many tools available for reading and writing JSON files
([Julia](https://github.com/JuliaLang/JSON.jl), [Python](https://docs.python.org/2.7/library/json.html)).
These tools can be used to generate a sky model using your favorite programming language.

**Fields:**

* `ref` is not used by TTCal, but is intended to be a record of the origin
  of the information used to define the source
* `name` is the name of the source
* `ra` and `dec` define the J2000 location of the source
* `I`, `Q`, `U`, and `V` define the Stokes parameters (in Jy) at the frequency `freq` (in Hz)
  `Q`, `U`, and `V` are optional and will be assumed to be zero if they are not given
* `index` defines the spectral index of the source (and higher order terms)

Higher order terms are defined such that

\[
    \log_{10} I = \log_{10} I_0 + \sum_{n=1}^N \alpha_n \log_{10}\left(\frac{\nu}{\nu_0}\right)^n
\]

where $I_0$ is the Stokes-I flux at the reference frequency and $\alpha_n$ represents
the values contained in `index`.

!!! note
    The Sun, the Moon, and Jupiter are special cases. These solar system objects do not
    need to have their position specified. TTCal will automatically determine their
    location based on the current time. This works based on the name of the source.
    If you name a source "Sun", "Moon", or "Jupiter" it will be placed at the
    correct location.

## Gaussian Sources

Gaussians sources can be added to the sky model by defining `major-fwhm`, `minor-fwhm`,
and `position-angle`.

* `major-fwhm` is the full width at half maximum (FWHM) of the major axis in arcseconds
* `minor-fwhm` is the FWHM of the minor axis in arcseconds
* `position-angle` is the angle the major axis makes with J2000 north where a positive angle is east of north

Note that the Stokes parameters `I`, `Q`, `U`, and `V` now give the integrated flux
over the entire source. An example Gaussian source is given below:

```
    ...
    {
        "ref": "Michael's imagination",
        "name": "A Gaussian source",
        "ra": "12h34m56.78s",
        "dec": "+12d34m56.78s",
        "I": 123.45,
        "freq": 70e6,
        "index": [-0.770],
        "major-fwhm": 200,
        "minor-fwhm": 100,
        "position-angle": -70
    }
    ...
```

## Multi-Component Sources

When peeling sources with TTCal each source receives its own peeling solution. This is not
always the desired behavior. For example if you are trying to peel Cyg A you might want to
use two sources to describe the two radio lobes, but both sources should be peeled with
the same solution. Multi-component sources allow you to define sources that are composed
of any number of the other source types. For example:

```
    ...
    {
        "name": "Cyg A",
        "components": [
            {
                "name": "east lobe",
                ...
            },
            {
                "name": "west lobe",
                ...
            }
        ]
    }
    ...
```

## RFI Sources

TTCal allows the possibility that some sources may not originate from the far-field of
the interferometer. That is, some sources may be so close to the interferometer that
the curvature of the incoming wavefront cannot be neglected. The exact transition between
the near-field and far-field depends on the wavelength, the maximum baseline of the
inteferometer, and the desired dynamic range. In order to correctly account for near-field
effects TTCal needs the longitude, latitude, and elevation of the source.

!!! note
    The ITRF coordinate system requires the distance to the center of the Earth in place of
    the elevation.

Sources that do live in the near-field are likely sources of radio frequency interference (RFI),
and therefore do not have the smooth power-law spectra expected of astronomical sources.
Specifying the flux of an RFI source therefore requires you to give a list of frequencies and
the flux at that frequency. For example:

```
    ...
    {
        "name": "Noisy Power Lines"
        "sys": "WGS84",
        "long": -118.31478,
        "lat": 37.14540,
        "el": 1226.709,
        "rfi-frequencies": [3.0e7, 3.75e7, 4.5e7, 5.25e7, 6.0e7],
        "rfi-I": [1.0, 2.0, 3.0, 4.0, 5.0]
    }
    ...
```

The longitude and latitude must be given in degrees and the elevation (or distance from
the center of the Earth) must be given in meters. The only valid coordinate systems are
currently `WGS84` and `ITRF`.
You may optionally also supply the fields `rfi-Q`, `rfi-U`, and `rfi-V` to specify the
corresponding Stokes parameters for the RFI emission.

