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
    getspec(visibilities, meta, point)
end

function getspec_internal!(output, visibilities, meta, source::Source)
    frame = reference_frame(meta)
    if isabovehorizon(frame, source)
        model = genvis_internal(meta, ConstantBeam(), [source])
        flatten_spectrum!(meta, model, source)
        getspec_internal!(output, visibilities, meta, model)
        # TEMP
        for β = 1:Nfreq(meta), α = 1:Nbase(meta)
            antenna1 = meta.baselines[α].antenna1
            antenna2 = meta.baselines[α].antenna2
            antenna1 == antenna2 && continue # don't use autocorrelations
            visibilities.flags[α,β] && continue
            V = visibilities.data[α,β]
            flux = get_total_flux(source, meta.channels[β])
        end
    end
end

function getspec_internal!(spectrum, visibilities, meta, model::Matrix{JonesMatrix})
    for β = 1:Nfreq(meta)
        numerator   = zero(JonesMatrix)
        denominator = zero(JonesMatrix)
        for α = 1:Nbase(meta)
            antenna1 = meta.baselines[α].antenna1
            antenna2 = meta.baselines[α].antenna2
            antenna1 == antenna2 && continue # don't use autocorrelations
            visibilities.flags[α,β] && continue
            numerator   += model[α,β]'*visibilities.data[α,β]
            denominator += model[α,β]'*model[α,β]

        end
        if abs(det(denominator)) > eps(Float64)
            spectrum[β] = make_hermitian(denominator \ numerator)
        end
    end
    spectrum
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

function flatten_spectrum!(meta, model, source::Source)
    for β = 1:Nfreq(meta)
        flux = get_total_flux(source, meta.channels[β])
        for α = 1:Nbase(meta)
            model[α,β] = model[α,β] / flux
            # Note that the following form is incorrect
            #     model[α,β] = flux \ model[α,β]
            # This decision is determined by the requirement that if the data and
            # model are exactly equal, we should get `flux` back exactly. The form
            # of the estimator in `getspec_internal!` (more specifically the order
            # of the matrix multiplications) then determines the order of the matrix
            # multiplications here.
        end
    end
end

function get_total_flux(source::Source, ν)
    linear(source.spectrum(ν))
end

function get_total_flux(source::MultiSource, ν)
    output = zero(HermitianJonesMatrix)
    for component in source.components
        output += get_total_flux(component, ν)
    end
    output
end

