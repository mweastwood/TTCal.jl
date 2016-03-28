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

Support for Gaussian sources is coming!

## Near-Field Sources

Experimental support for near-field sources is currently available in the latest version
of TTCal. Because this is an experimental feature it is not yet very easy to use.

## Diffuse Emission

Diffuse emission is currently not supported within TTCal. Eventual support is on the
horizon.

