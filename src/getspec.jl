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

"""
    getspec(visibilities, metadata, source)
    getspec(visibilities, metadata, direction)

*Description*

This function will measure the flux as a function of frequency for the given source
or the given direction. All four correlations are included in this measurement so
that all four Stokes parameters can be calculated if desired.

The flux of a resolved source will be systematically underestimated by this routine
if you only supply a direction. However if you provide the full source model
including the resolved components then `getspec` will correctly account for the
resolved flux in its measurement.

Note that this function is susceptible to sidelobe contamination from other
bright sources in the field.

*Arguments*

* `visibilities` - the measured visibilities from which to estimate the flux
* `metadata` - the metadata describing the interferometer
* `source` or `direction` - the source model or direction in which to estimate the flux

Note that if you supply a direction, `getspec` will estimate the flux of a point
source in that direction. If you know your source is resolved, supply a source model
instead.
"""
function getspec(visibilities::Visibilities, meta::Metadata, source::Source)
    spectrum = zeros(HermitianJonesMatrix, Nfreq(meta))
    getspec_internal!(spectrum, visibilities, meta, source)
    spectrum
end

function getspec(visibilities::Visibilities, meta::Metadata, direction::Direction)
    flat = PowerLaw(1, 0, 0, 0, 10e6, [0.0])
    point = PointSource("dummy", direction, flat)
    getspec(visibilities, meta, source)
end

function getspec_internal!(spectrum, visibilities, meta::Metadata, source::Source)
    # TODO account for the spectrum of the source
    if isabovehorizon(frame, source)
        model = genvis_internal(meta, ConstantBeam(), source)
        getspec_internal!(spectrum, visibilities, model)
    end
end

function getspec_internal!(spectrum, visibilities, model::Matrix{JonesMatrix})
    for β = 1:Nfreq(visibilities)
    end
end

#function getspec_internal!(spectrum, visibilities, meta, direction, phase_center)
#    delays = geometric_delays(meta.antennas, direction, phase_center)
#    for β = 1:Nfreq(meta)
#        count = 0 # the number of baselines used in the calculation
#        frequency = meta.channels[β]
#        fringes = delays_to_fringes(delays, frequency)
#        for α = 1:Nbase(meta)
#            antenna1 = meta.baselines[α].antenna1
#            antenna2 = meta.baselines[α].antenna2
#            antenna1 == antenna2 && continue # don't use autocorrelations
#            visibilities.flags[α,β] && continue
#            fringe = conj(fringes[antenna1]) * fringes[antenna2]
#            xx = real(visibilities.data[α,β].xx*fringe)
#            xy = complex(real(visibilities.data[α,β].xy*fringe), -imag(visibilities.data[α,β].yx*fringe))
#            yy = real(visibilities.data[α,β].yy*fringe)
#            spectrum[β] += HermitianJonesMatrix(xx, xy, yy)
#            count += 1
#        end
#        if count > 0
#            spectrum[β] /= count
#        else
#            spectrum[β] = HermitianJonesMatrix(0, 0, 0)
#        end
#    end
#    spectrum
#end

#function with_flat_spectrum(source::Source)
#    mysource = deepcopy(source) # don't modify the argument!
#    mysource.spectrum
#end

