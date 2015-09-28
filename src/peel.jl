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
                          ms::Table,
                          sources::Vector{PointSource};
                          maxiter = 20, tolerance = 1e-3,
                          minuvw = 0.0)

Peel the given list of sources from the measurement set.

The type supplied as the first argument determines the
manner in which the sources are peeled:

* `PolarizationCalibration` - each source receives a full set of Jones matrices
* `GainCalibration` - each source receives a full set of complex gains
* `AmplitudeCalibration` - each source receives a full set of gain amplitudes
"""
function peel!{T<:Calibration}(::Type{T},
                               ms::Table,
                               sources::Vector{PointSource};
                               maxiter::Int = 20,
                               tolerance::Float64 = 1e-3,
                               minuvw::Float64 = 0.0)
    phase_dir = MeasurementSets.phase_direction(ms)
    u,v,w = MeasurementSets.uvw(ms)
    ν = MeasurementSets.frequency(ms)
    ant1,ant2 = MeasurementSets.antennas(ms)

    frame = ReferenceFrame()
    set!(frame,MeasurementSets.position(ms))
    set!(frame,MeasurementSets.time(ms))
    sources = filter(source -> isabovehorizon(frame,source),sources)

    Nant    = numrows(Table(ms[kw"ANTENNA"]))
    Nfreq   = length(ν)

    data  = MeasurementSets.corrected_data(ms)
    flags = MeasurementSets.flags(ms)

    calibrations = [T(Nant,Nfreq) for source in sources]
    coherencies  = [genvis(frame,phase_dir,[source],u,v,w,ν) for source in sources]

    peel!(calibrations,coherencies,data,flags,
          u,v,w,ν,ant1,ant2,maxiter,tolerance,minuvw)
    ms["CORRECTED_DATA"] = data
    calibrations
end

function peel!(calibrations,coherencies,data,flags,
               u,v,w,ν,ant1,ant2,maxiter,tolerance,minuvw)
    Nsource = length(calibrations)
    Nfreq = size(data,2)
    Nbase = size(data,3)

    # Subtract all of the sources
    # (assuming the beam is unity towards each source)
    for coherency in coherencies
        subsrc!(data,coherency)
    end

    # Flag all of the short baselines
    # (so that they do not contribute to the fit)
    for α = 1:Nbase, β = 1:Nfreq
        if sqrt(u[α]^2 + v[α]^2 + w[α]^2) < minuvw*c/ν[β]
            flags[:,β,α] = true
        end
    end

    # Derive a calibration towards each source
    for iter = 1:3
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
                   ant1,ant2,maxiter,tolerance)

            # Take the source back out of the measured visibilities,
            # but this time subtract it with the corrected gains toward
            # that direction.
            corrupted = copy(coherency)
            corrupt!(corrupted,calibration_toward_source,ant1,ant2)
            subsrc!(data,corrupted)
        end
    end
    nothing
end

