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
    peel!{T<:Calibration}(::Type{T}, visibilities::Visibilities, meta::Metadata, sources;
                          peeliter = 3, maxiter = 20, tolerance = 1e-3)

Peel the given list of sources from the measurement set.

The type supplied as the first argument determines the
manner in which the sources are peeled:

* `PolarizationCalibration` - each source receives a full set of Jones matrices
* `GainCalibration` - each source receives a full set of complex gains
* `AmplitudeCalibration` - each source receives a full set of gain amplitudes
"""
function peel!{T<:Calibration}(::Type{T}, visibilities::Visibilities, meta::Metadata, sources;
                               peeliter::Int = 3, maxiter::Int = 20, tolerance::Float64 = 1e-3)
    frame = reference_frame(meta)
    sources = abovehorizon(frame, sources)
    calibrations = [T(Nant(meta), Nfreq(meta)) for source in sources]
    coherencies  = [genvis(meta, source) for source in sources]
    peel!(calibrations, coherencies, visibilities, meta, peeliter, maxiter, tolerance)
    calibrations
end

function peel!(calibrations, coherencies, visibilities, meta, peeliter, maxiter, tolerance)
    Nsource = length(calibrations)

    # Subtract all of the sources
    # (assuming the beam is unity towards each source)
    for coherency in coherencies
        subsrc!(visibilities, coherency)
    end

    # Derive a calibration towards each source
    p = Progress(peeliter*Nsource)
    for iter = 1:peeliter
        for s = 1:Nsource
            coherency = coherencies[s]
            calibration_toward_source = calibrations[s]

            # Put one source back into the visibilities.
            corrupted = deepcopy(coherency)
            corrupt!(corrupted, meta, calibration_toward_source)
            putsrc!(visibilities, corrupted)

            # Solve for the calibration in the direction
            # of the current source.
            solve!(calibration_toward_source, visibilities, coherency,
                   meta, maxiter, tolerance, quiet = true)

            # Take the source back out of the measured visibilities,
            # but this time subtract it with the corrected gains toward
            # that direction.
            corrupted = deepcopy(coherency)
            corrupt!(corrupted, meta, calibration_toward_source)
            subsrc!(visibilities, corrupted)

            next!(p)
        end
    end
    for calibration in calibrations
        if sum(calibration.flags) > 0.5length(calibration.flags)
            warn("Frequently failed to converge. There will likely be large residuals.")
        end
    end
    calibrations
end

