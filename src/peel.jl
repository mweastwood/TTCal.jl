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
    peel!(visibilities, metadata, sources)

*Description*

Peel sources from the visibilities.

*Arguments*

* `visibilities` - a calibrated set of visibilities
* `metadata` - the metadata describing the interferometer
* `sources` - the list of sources to peel from the visibilities

*Keyword Arguments*

* `peeliter` - the number of peeling iterations
* `maxiter` - the maximum number of iterations to take on each frequency channel (defaults to `20`)
* `tolerance` - the relative tolerance used to test for convergence (defaults to `1e-3`)
* `quiet` - suppresses printing if set to `true` (defaults to `false`)

Note that peeling iterates through the list of sources in order. A peeling
iteration is defined as one pass through the full list. Usually 2
or 3 iterations seems to be sufficient.
"""
function peel!(visibilities::Visibilities, meta::Metadata, sources;
               peeliter = 3, maxiter = 20, tolerance = 1e-3, quiet = false)
    frame = reference_frame(meta)
    sources = abovehorizon(frame, sources)
    calibrations = [GainCalibration(Nant(meta), Nfreq(meta)) for source in sources]
    coherencies  = [genvis(meta, source) for source in sources]
    peel!(calibrations, coherencies, visibilities, meta, peeliter, maxiter, tolerance, quiet)
    calibrations
end

"""
    shave!(visibilities, metadata, sources)

*Description*

Shave sourcces from the visibilities. This is essentially equivalent
to peeling, but only one calibration is applied across the entire
frequency channel (instead of one calibration per frequency channel).

Note that the term "shave" is non-standard. I picked the name because
it's short, a synonym for "peel", and will probably force you to read
the documentation to figure out how `shave!` differs from `peel!`.
Looks like it worked!

*Arguments*

* `visibilities` - a calibrated set of visibilities
* `metadata` - the metadata describing the interferometer
* `sources` - the list of sources to peel from the visibilities

*Keyword Arguments*

* `peeliter` - the number of peeling iterations
* `maxiter` - the maximum number of iterations to take on each frequency channel (defaults to `20`)
* `tolerance` - the relative tolerance used to test for convergence (defaults to `1e-3`)
* `quiet` - suppresses printing if set to `true` (defaults to `false`)

Note that peeling iterates through the list of sources. A peeling
iteration is defined as one pass through the full list. Usually 2
or 3 iterations seems to be sufficient.
"""
function shave!(visibilities::Visibilities, meta::Metadata, sources;
                peeliter = 3, maxiter = 20, tolerance = 1e-3, quiet = false)
    frame = reference_frame(meta)
    sources = abovehorizon(frame, sources)
    calibrations = [GainCalibration(Nant(meta), 1) for source in sources]
    coherencies  = [genvis(meta, source) for source in sources]
    peel!(calibrations, coherencies, visibilities, meta, peeliter, maxiter, tolerance, quiet)
    calibrations
end

function peel!(calibrations, coherencies, visibilities, meta, peeliter, maxiter, tolerance, quiet)
    Nsource = length(calibrations)

    # Subtract all of the sources
    # (assuming the beam is unity towards each source)
    for coherency in coherencies
        subsrc!(visibilities, coherency)
    end

    # Derive a calibration towards each source
    quiet || (p = Progress(peeliter*Nsource, "Peeling: "))
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
            do_peeling_solve!(calibration_toward_source, visibilities, coherency, meta, maxiter, tolerance)

            # Take the source back out of the measured visibilities,
            # but this time subtract it with the corrected gains toward
            # that direction.
            corrupted = deepcopy(coherency)
            corrupt!(corrupted, meta, calibration_toward_source)
            subsrc!(visibilities, corrupted)

            quiet || next!(p)
        end
    end
    if !quiet
        for calibration in calibrations
            if sum(calibration.flags) > 0.5length(calibration.flags)
                warn("Frequently failed to converge. There will likely be large residuals.")
            end
        end
    end
    calibrations
end

function do_peeling_solve!(calibration_toward_source, visibilities, coherency, meta, maxiter, tolerance)
    if Nfreq(calibration_toward_source) == 1
        # shaving
        solve_allchannels!(calibration_toward_source, visibilities, coherency, meta, maxiter, tolerance)
    else
        # peeling
        solve!(calibration_toward_source, visibilities, coherency, meta, maxiter, tolerance, true)
    end
end

