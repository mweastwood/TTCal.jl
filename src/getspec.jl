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

#"""
#    getspec(visibilities, metadata, source)
#    getspec(visibilities, metadata, direction)
#
#*Description*
#
#This function will measure the flux as a function of frequency for the given source
#or the given direction. All four correlations are included in this measurement so
#that all four Stokes parameters can be calculated if desired.
#
#The flux of a resolved source will be systematically underestimated by this routine
#if you only supply a direction. However if you provide the full source model
#including the resolved components then `getspec` will correctly account for the
#resolved flux in its measurement.
#
#Note that this function is susceptible to sidelobe contamination from other
#bright sources in the field.
#
#*Arguments*
#
#* `visibilities` - the measured visibilities from which to estimate the flux
#* `metadata` - the metadata describing the interferometer
#* `source` or `direction` - the source model or direction in which to estimate the flux
#
#Note that if you supply a direction, `getspec` will estimate the flux of a point
#source in that direction. If you know your source is resolved, supply a source model
#instead.
#"""
function getspec(dataset::Dataset, source::Union{Source, Direction, Position})
    spectrum = zeros(StokesVector, Nfreq(dataset))
    getspec_internal!(spectrum, dataset, source)
    spectrum
end

function getflux(dataset::Dataset, source::Union{Source, Direction, Position})
    mean(getspec(dataset, source))
end

function getspec_internal!(output, dataset, source::Source)
    # TODO: handle multiple time integrations correctly
    frame = ReferenceFrame(dataset.metadata)
    if isabovehorizon(frame, source)
        model = genvis(dataset.metadata, ConstantBeam(), source, polarization=polarization(dataset))
        flatten_spectrum!(model, source)
        getspec_internal!(output, dataset, model)
    end
end

function getspec_internal!(spectrum, dataset, model::Dataset)
    for frequency = 1:Nfreq(dataset)
        visibilities       = dataset[frequency, 1]
        model_visibilities =   model[frequency, 1]
        numerator   = zero(eltype(dataset))
        denominator = zero(eltype(dataset))
        for ant1 = 1:Nant(dataset), ant2 = ant1+1:Nant(dataset)
            isflagged(visibilities, ant1, ant2) && continue
            J1 = visibilities[ant1, ant2]
            J2 = model_visibilities[ant1, ant2]
            J′ = J2'
            numerator   += J′*J1
            denominator += J′*J2
        end
        if abs(det(denominator)) > eps(Float64)
            spectrum[frequency] = create_stokes_vector(denominator \ numerator,
                                                       polarization(dataset))
        end
    end
    spectrum
end

function getspec_internal!(spectrum, dataset, direction::Direction)
    # TODO: handle multiple time integrations correctly
    metadata = dataset.metadata
    frame = ReferenceFrame(metadata)
    itrf_direction    = measure(frame, direction, dir"ITRF")
    itrf_phase_center = measure(frame, metadata.phase_centers[1], dir"ITRF")
    delays = geometric_delays(metadata.positions, itrf_direction, itrf_phase_center)
    getspec_from_delays!(spectrum, dataset, delays)
end

function getspec_internal!(spectrum, dataset, position::Position)
    # TODO: handle multiple time integrations correctly
    metadata = dataset.metadata
    frame = ReferenceFrame(metadata)
    itrf_position     = measure(frame, position, pos"ITRF")
    itrf_phase_center = measure(frame, metadata.phase_centers[1], dir"ITRF")
    delays = geometric_delays(metadata.positions, itrf_position, itrf_phase_center)
    getspec_from_delays!(spectrum, dataset, delays)
end

function getspec_from_delays!(spectrum, dataset, delays)
    metadata = dataset.metadata
    for frequency = 1:Nfreq(dataset)
        fringes = delays_to_fringes(delays, metadata.frequencies[frequency])
        visibilities = dataset[frequency, 1]
        numerator   = zero(eltype(dataset))
        denominator = 0
        for antenna1 = 1:Nant(visibilities), antenna2 = antenna1+1:Nant(visibilities)
            if !isflagged(visibilities, antenna1, antenna2)
                fringe = fringes[antenna1] * conj(fringes[antenna2])
                numerator   += conj(fringe)*visibilities[antenna1, antenna2]
                denominator += 1
            end
        end
        if denominator > 0
            spectrum[frequency] = create_stokes_vector(denominator \ numerator,
                                                       polarization(dataset))
        end
    end
    spectrum
end

create_stokes_vector(value, ::Any) = StokesVector(make_hermitian(value))
function create_stokes_vector(value, ::Type{TTCal.XX})
    StokesVector(HermitianJonesMatrix(real(value), 0, 0))
end
function create_stokes_vector(value, ::Type{TTCal.YY})
    StokesVector(HermitianJonesMatrix(0, 0, real(value)))
end

# Design justification
#
# Let's say we are trying to measure the flux of a multi-component source.
# The user has likely gone to great pains to get the relative fluxes of the
# various components just right. We do not want to discard the spectral
# information encoded in all these various components.
#
# So we generate the model visibilities using the input source model as given.
# We then scale these visibilities by the total flux of the source (counting
# all of the components) so that the flux of the model visibilities is unity.

