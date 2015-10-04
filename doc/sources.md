## Sky Models

TTCal expects the sky model to be specified as a list of point sources in a JSON file.
Each source receives a name (`name`), J2000 position (`ra` and `dec`), and power-law spectrum.
The Stokes parameters (`I`, `Q`, `U`, and `V`) must be given at the reference frequency
(`freq`, measured in Hz). The spectral index is specified by `index`. Higher order terms
may also be specified. See the following example for Cas A and Cyg A using the published
spectra for Baars et al. 1977.

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

