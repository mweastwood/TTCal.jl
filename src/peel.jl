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
    peel!{T<:Calibration}(::Type{T},
                          ms::MeasurementSet,
                          sources::Vector{Source},
                          beam::BeamModel;
                          peeliter = 3,
                          maxiter = 20,
                          tolerance = 1e-3,
                          minuvw = 0.0)

Peel the given list of sources from the measurement set.

The type supplied as the first argument determines the
manner in which the sources are peeled:

* `PolarizationCalibration` - each source receives a full set of Jones matrices
* `GainCalibration` - each source receives a full set of complex gains
* `AmplitudeCalibration` - each source receives a full set of gain amplitudes
"""
function peel!{T<:Calibration}(::Type{T},
                               ms::MeasurementSet,
                               sources::Vector{Source},
                               beam::BeamModel;
                               peeliter::Int = 3,
                               maxiter::Int = 20,
                               tolerance::Float64 = 1e-3,
                               minuvw::Float64 = 0.0)
    sources = abovehorizon(ms.frame,sources)
    data  = get_corrected_data(ms)
    flags = get_flags(ms)
    flag_short_baselines!(flags,minuvw,ms.u,ms.v,ms.w,ms.ν)
    calibrations = [T(ms.Nant,ms.Nfreq) for source in sources]
    coherencies  = [genvis(ms,source,beam) for source in sources]
    peel!(calibrations,coherencies,data,flags,
          ms.u,ms.v,ms.w,ms.ν,ms.ant1,ms.ant2,
          peeliter,maxiter,tolerance)
    set_corrected_data!(ms,data,true)
    calibrations
end

function peel!(calibrations,coherencies,data,flags,
               u,v,w,ν,ant1,ant2,peeliter,maxiter,tolerance)
    Nsource = length(calibrations)
    Nfreq = size(data,2)
    Nbase = size(data,3)

    # Subtract all of the sources
    # (assuming the beam is unity towards each source)
    for coherency in coherencies
        subsrc!(data,coherency)
    end

    # Derive a calibration towards each source
    p = Progress(peeliter*Nsource, 1, "Peeling...", 50)
    for iter = 1:peeliter
        for s = 1:Nsource
            coherency = coherencies[s]
            calibration_toward_source = calibrations[s]

            # Put one source back into the visibilities.
            corrupted = copy(coherency)
            corrupt!(corrupted,calibration_toward_source,ant1,ant2)
            putsrc!(data,corrupted)

            # Solve for the calibration in the direction
            # of the current source.
            solve!(calibration_toward_source,
                   data,coherency,flags,
                   ant1,ant2,maxiter,tolerance,
                   quiet = true)

            # Take the source back out of the measured visibilities,
            # but this time subtract it with the corrected gains toward
            # that direction.
            corrupted = copy(coherency)
            corrupt!(corrupted,calibration_toward_source,ant1,ant2)
            subsrc!(data,corrupted)

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

