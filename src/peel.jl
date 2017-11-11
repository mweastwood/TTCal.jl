# Copyright (c) 2015, 2016 Michael Eastwood
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

"A wrapper around `Source` that tells peeling which algorithm to use."
abstract AbstractPeelingSource

macro peelingsource(name, calibration_type_expr)
    expr = quote
        Base.@__doc__ immutable $name <: AbstractPeelingSource
            source :: Source
        end
        calibration_type(meta::Metadata, ::$name) = $calibration_type_expr
    end
    esc(expr)
end

@inline unwrap(source::AbstractPeelingSource) = source.source
@inline unwrap(source) = source

"One diagonal Jones matrix per antenna per frequency channel."
@peelingsource PeelingSource GainCalibration(Nant(meta), Nfreq(meta))

"One diagonal Jones matrix per antenna per subband."
@peelingsource ShavingSource GainCalibration(Nant(meta), 1)

"One full Jones matrix per antenna per frequency channel."
@peelingsource ZestingSource PolarizationCalibration(Nant(meta), Nfreq(meta))

"One full Jones matrix per antenna per subband."
@peelingsource PruningSource PolarizationCalibration(Nant(meta), 1)

"""
    peel!(visibilities, metadata, beam, sources)

*Description*

Peel sources from the visibilities.

| Name    | Calibration Type      | Bandwidth       |
|---------|-----------------------|-----------------|
| Peeling | Diagonal Jones matrix | One per channel |
| Shaving | Diagonal Jones matrix | One per subband |
| Zesting | Full Jones matrix     | One per channel |
| Pruning | Full Jones Matrix     | One per subband |

Note that the names "shaving", "zesting", and "pruning" are non-standard
puns. I picked these names because they are short synonyms for "peeling"
and will likely force you to check the documentation when you come across
them. Looks like it worked!

*Arguments*

* `visibilities` - a calibrated set of visibilities
* `metadata` - the metadata describing the interferometer
* `beam` - the primary beam model
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
function peel!{T<:AbstractPeelingSource}(visibilities::Visibilities, meta::Metadata,
                                         beam::BeamModel, sources::Vector{T};
                                         peeliter = 3, maxiter = 20, tolerance = 1e-3, quiet = false)
    frame = reference_frame(meta)
    sources = filter(sources) do source
        isabovehorizon(frame, unwrap(source))
    end
    calibrations = [calibration_type(meta, source) for source in sources]
    coherencies  = [genvis(meta, beam, unwrap(source)) for source in sources]
    peel!(calibrations, coherencies, visibilities, meta, peeliter, maxiter, tolerance, quiet)
    calibrations
end

function peel!(calibrations, coherencies, visibilities, meta, peeliter, maxiter, tolerance, quiet)
    Nsource = length(calibrations)

    # Subtract all of the sources
    # (assuming the beam is unity towards each source)
    for s = 1:Nsource
        # Make sure to apply the initial gains here because we're going to be restoring the source
        # with those same initial gains
        coherency = coherencies[s]
        calibration_toward_source = calibrations[s]
        corrupted = deepcopy(coherency)
        corrupt!(corrupted, meta, calibration_toward_source)
        subsrc!(visibilities, corrupted)
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
        for (idx, calibration) in enumerate(calibrations)
            if sum(calibration.flags) > 0.5length(calibration.flags)
                warn("Calibration number $idx frequently failed to converge.")
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

# Legacy methods

function peel!{T<:Source}(visibilities::Visibilities, metadata::Metadata, beam::BeamModel, sources::Vector{T};
                          peeliter = 3, maxiter = 20, tolerance = 1e-3, quiet = false)
    peelingsources = [PeelingSource(source) for source in sources]
    peel!(visibilities, metadata, beam, peelingsources,
          peeliter=peeliter, maxiter=maxiter, tolerance=tolerance, quiet=quiet)
end

function shave!{T<:Source}(visibilities::Visibilities, metadata::Metadata, beam::BeamModel, sources::Vector{T};
                           peeliter = 3, maxiter = 20, tolerance = 1e-3, quiet = false)
    shavingsources = [ShavingSource(source) for source in sources]
    peel!(visibilities, metadata, beam, shavingsources,
          peeliter=peeliter, maxiter=maxiter, tolerance=tolerance, quiet=quiet)
end

