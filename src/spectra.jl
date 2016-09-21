# Copyright (c) 2015 Michael Eastwood
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

abstract Spectrum

doc"""
    PowerLaw <: Spectrum

A multi-component power-law spectrum.

``` math
    \log_{10} S = \log_{10} S_0 + \sum_{n=1}^N \alpha_n \log_{10}\left(\frac{\nu}{\nu_0}\right)^n
```

where $S$ is a Stokes parameter, and $\alpha_n$ is the list of spectral indices. At least one
spectral index needs to be provided.
"""
type PowerLaw <: Spectrum
    stokes :: StokesVector
    ν :: Float64 # Hz
    α :: Vector{Float64}
end

PowerLaw(I, Q, U, V, ν, index) = PowerLaw(StokesVector(I, Q, U, V), ν, index)
PowerLaw(I::Real, ν, index) = PowerLaw(StokesVector(I, 0, 0, 0), ν, index)

function (spectrum::PowerLaw)(ν::AbstractFloat)
    s = sign(spectrum.stokes.I)
    log_I = log10(abs(spectrum.stokes.I))
    log_ν = log10(ν/spectrum.ν)
    for (i,α) in enumerate(spectrum.α)
        log_I += α*log_ν^i
    end
    I = s*10^log_I
    Q = spectrum.stokes.Q / spectrum.stokes.I * I
    U = spectrum.stokes.U / spectrum.stokes.I * I
    V = spectrum.stokes.V / spectrum.stokes.I * I
    StokesVector(I, Q, U, V)
end

function ==(lhs::PowerLaw, rhs::PowerLaw)
    lhs.stokes == rhs.stokes && lhs.ν == rhs.ν && lhs.α == rhs.α
end

"""
    RFISpectrum <: Spectrum

A simple list of frequency channels and the Stokes parameters at each channel.
"""
type RFISpectrum <: Spectrum
    channels :: Vector{Float64}
    stokes   :: Vector{StokesVector}
    function RFISpectrum(channels, stokes)
        if length(channels) != length(stokes)
            error("number of frequency channels must match number of Stokes vectors")
        end
        new(channels, stokes)
    end
end

function Base.show(io::IO, spectrum::RFISpectrum)
    N = length(spectrum.channels)
    ν1 = spectrum.channels[1] / 1e6
    ν2 = spectrum.channels[2] / 1e6
    stokes = mean(spectrum.stokes)
    print(io, @sprintf("RFISpectrum(%d channels between %.3f and %.3f MHz", N, ν1, ν2), ", ", stokes, ")")
end

function (spectrum::RFISpectrum)(ν::AbstractFloat)
    idx = indmin(abs2(channel - ν) for channel in spectrum.channels)
    spectrum.stokes[idx]
end

function ==(lhs::RFISpectrum, rhs::RFISpectrum)
    lhs.channels == rhs.channels && lhs.stokes == rhs.stokes
end

