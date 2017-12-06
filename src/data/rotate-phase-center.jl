# Copyright (c) 2015-2017 Michael Eastwood
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

function rotate_phase_center!(dataset::Dataset, new_phase_center::Direction)
    metadata = dataset.metadata
    frame = ReferenceFrame(metadata)
    for (idx, time) in enumerate(metadata.times)
        set!(frame, time)
        old_phase_center_itrf = measure(frame, metadata.phase_centers[idx], dir"ITRF")
        new_phase_center_itrf = measure(frame, new_phase_center, dir"ITRF")
        for (jdx, frequency) in enumerate(metadata.frequencies)
            visibilities = dataset[jdx, idx]
            rotate_phase_center_onechannel!(visibilities, frequency, metadata.positions,
                                            old_phase_center_itrf, new_phase_center_itrf)
        end
    end
    dataset
end

function rotate_phase_center_onechannel!(visibilities, frequency, positions,
                                         old_phase_center, new_phase_center)
    delays = geometric_delays(positions, old_phase_center, new_phase_center)
    fringes = delays_to_fringes(delays, frequency)
    for antenna1 = 1:Nant(visibilities), antenna2 = antenna1:Nant(visibilities)
        fringe = fringes[antenna1] * conj(fringes[antenna2])
        visibilities[antenna1, antenna2] *= fringe
    end
    visibilities
end

