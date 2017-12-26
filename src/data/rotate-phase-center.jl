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
        dataset.metadata.phase_centers[idx] = new_phase_center
    end
    dataset
end

function rotate_phase_center_onechannel!(visibilities, frequency, positions,
                                         old_phase_center, new_phase_center)
    delays = geometric_delays(positions, old_phase_center, new_phase_center)
    fringes = delays_to_fringes(delays, frequency)
    for antenna1 = 1:Nant(visibilities), antenna2 = antenna1:Nant(visibilities)
        if !isflagged(visibilities, antenna1, antenna2)
            fringe = fringes[antenna1] * conj(fringes[antenna2])
            visibilities[antenna1, antenna2] *= fringe
        end
    end
    visibilities
end

"Bring the source into focus."
function rotate_phase_center!(dataset::Dataset, source::Source)
    frame = ReferenceFrame(dataset.metadata)
    model = genvis(dataset.metadata, ConstantBeam(), source, polarization=polarization(dataset))
    flatten_spectrum!(model, source)
    for time = 1:Ntime(dataset)
        set!(frame, dataset.metadata.times[time])
        for frequency = 1:Nfreq(dataset)
            visibilities       = dataset[frequency, time]
            model_visibilities =   model[frequency, time]
            rotate_phase_center_onechannel!(visibilities, model_visibilities)
        end
        direction = mean_direction(frame, source, mean(dataset.metadata.frequencies))
        dataset.metadata.phase_centers[time] = direction
    end
end

function rotate_phase_center_onechannel!(data::Visibilities, model::Visibilities)
    for ant1 = 1:Nant(data), ant2 = ant1:Nant(data)
        # TODO is this order, correct? it doesn't matter if these are diagonal matrices...
        if !isflagged(data, ant1, ant2)
            data[ant1, ant2] = data[ant1, ant2] / model[ant1, ant2]
        end
    end
end

function flatten_spectrum!(model, source::Source)
    # TODO: handle multiple time integrations correctly
    for frequency = 1:Nfreq(model)
        visibilities = model[frequency, 1]
        flux = total_flux(source, model.metadata.frequencies[frequency])
        jones = HermitianJonesMatrix(flux)
        for ant1 = 1:Nant(model), ant2 = ant1:Nant(model)
            # Note that the following form is incorrect
            #     model[ant1, ant2] = flux \ model[ant1, ant2]
            # This decision is determined by the requirement that if the data and
            # model are exactly equal, we should get `flux` back exactly. The form
            # of the estimator in `getspec_internal!` (more specifically the order
            # of the matrix multiplications) then determines the order of the matrix
            # multiplications here.
            if !isflagged(visibilities, ant1, ant2)
                visibilities[ant1, ant2] = flatten_spectrum(visibilities[ant1, ant2], jones,
                                                            polarization(model))
            end
        end
    end
end

flatten_spectrum(value, jones, ::Any) = value/jones
flatten_spectrum(value, jones, ::Type{TTCal.XX}) = value/jones.xx
flatten_spectrum(value, jones, ::Type{TTCal.YY}) = value/jones.yy

