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
    getspec(visibilities::Visibilities, meta::Metadata, direction::Direction)

This function measures the flux as a function of frequency in a given direction.

Note that this function is susceptible to sidelobe contamination from other
bright sources in the field.
"""
function getspec(visibilities::Visibilities, meta::Metadata, direction::Direction)
    frame = reference_frame(meta)
    spectrum = zeros(HermitianJonesMatrix, Nfreq(meta))
    if isabovehorizon(frame, direction)
        direction = measure(frame, direction, dir"ITRF")
        phase_center = measure(frame, meta.phase_center, dir"ITRF")
        getspec_internal!(spectrum, visibilities, meta, direction, phase_center)
    end
    spectrum
end

function getspec_internal!(spectrum, visibilities, meta, direction, phase_center)
    delays = geometric_delays(meta.antennas, direction, phase_center)
    for β = 1:Nfreq(meta)
        count = 0 # the number of baselines used in the calculation
        frequency = meta.channels[β]
        fringes = delays_to_fringes(delays, frequency)
        for α = 1:Nbase(meta)
            antenna1 = meta.baselines[α].antenna1
            antenna2 = meta.baselines[α].antenna2
            antenna1 == antenna2 && continue # don't use autocorrelations
            visibilities.flags[α,β] && continue
            fringe = conj(fringes[antenna1]) * fringes[antenna2]
            xx = real(visibilities.data[α,β].xx*fringe)
            xy = complex(real(visibilities.data[α,β].xy*fringe), -imag(visibilities.data[α,β].yx*fringe))
            yy = real(visibilities.data[α,β].yy*fringe)
            spectrum[β] += HermitianJonesMatrix(xx, xy, yy)
            count += 1
        end
        spectrum[β] /= count
    end
    spectrum
end

