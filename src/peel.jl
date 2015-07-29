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

################################################################################
# Public Interface

function peel!(ms::Table,
               sources::Vector{PointSource};
               maxiter::Int = 10,
               tolerance::Float64 = 1e-3,
               minuvw::Float64 = 0.0)
    dir   = MeasurementSets.phase_direction(ms)
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

    calibrations = [GainCalibration(Nant,Nfreq) for source in sources]
    coherencies  = [genvis(dir,source,u,v,w,ν)  for source in sources]

    peel_internal!(calibrations,coherencies,data,flags,
                   u,v,w,ν,ant1,ant2,maxiter,tolerance,minuvw)
    ms["CORRECTED_DATA"] = data
    calibrations
end

################################################################################
# Internal Interface

function peel_internal!(calibrations,coherencies,data,flags,
                        u,v,w,ν,ant1,ant2,maxiter,tolerance,minuvw)
    Nsource = length(calibrations)
    Nfreq = size(data,2)
    Nbase = size(data,3)

    # Subtract all of the sources
    # (assuming the beam is unity towards each source)
    for coherency in coherencies
        for i in eachindex(data)
            data[i] -= coherency[i]
        end
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
            for i in eachindex(data)
                data[i] += corrupted[i]
            end

            bandpass_internal!(calibration_toward_source,
                               data,coherency,flags,
                               ant1,ant2,maxiter,tolerance,1)

            # Take the source back out of the measured visibilities,
            # but this time subtract it with the corrected gains toward
            # that direction.
            corrupted = copy(coherency)
            corrupt!(corrupted,calibration_toward_source,ant1,ant2)
            for i in eachindex(data)
                data[i] -= corrupted[i]
            end
        end
    end
    nothing
end

function corrupt!(data::Array{Complex64,3},
                  cal::GainCalibration,
                  ant1,ant2)
    Nbase = length(ant1)
    for α = 1:Nbase, β = 1:Nfreq(cal)
        data[1,β,α] = cal.gains[ant1[α],β,1]*conj(cal.gains[ant2[α],β,1])*data[1,β,α]
        data[2,β,α] = cal.gains[ant1[α],β,1]*conj(cal.gains[ant2[α],β,2])*data[2,β,α]
        data[3,β,α] = cal.gains[ant1[α],β,2]*conj(cal.gains[ant2[α],β,1])*data[3,β,α]
        data[4,β,α] = cal.gains[ant1[α],β,2]*conj(cal.gains[ant2[α],β,2])*data[4,β,α]
    end
end

