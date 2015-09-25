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
    AmplitudeCalibration <: Calibration

This type stores the information for calibrating only
the amplitude of the electronic gains. It stores
the gain amplitudes and flags for each antenna, frequency channel,
and polarization.
"""
immutable AmplitudeCalibration <: Calibration
    amplitudes::Array{Float64,3}
    flags::Array{Bool,3}
end

"""
    AmplitudeCalibration(Nant,Nfreq)

Create a calibration table for `Nant` antennas with
`Nfreq` frequency channels where all the amplitudes
are initially set to unity.
"""
function AmplitudeCalibration(Nant,Nfreq)
    amplitudes = ones(Float64,Nant,Nfreq,2)
    flags = zeros(Bool,Nant,Nfreq,2)
    AmplitudeCalibration(amplitudes,flags)
end

Nant(g::AmplitudeCalibration) = size(g.amplitudes,1)
Nfreq(g::AmplitudeCalibration) = size(g.amplitudes,2)

doc"""
    invert(cal::AmplitudeCalibration)

Returns the inverse of the given calibration.
The gain amplitude $a$ of each antenna is set to $1/a$.
"""
function invert(cal::AmplitudeCalibration)
    output = AmplitudeCalibration(Nant(cal),Nfreq(cal))
    for i in eachindex(cal.amplitudes)
        output.amplitudes[i] = inv(cal.amplitudes[i])
    end
    output
end

"""
    ampcal(ms::Table,sources::Vector{PointSource};
           maxiter = 30, tolerance = 1e-3,
           force_imaging_columns = false)

Solve for the amplitude of the interferometer's gains.
"""
function ampcal(ms::Table,
                sources::Vector{PointSource};
                maxiter::Int = 30,
                tolerance::Float64 = 1e-3,
                force_imaging_columns::Bool = false)
    phase_dir = MeasurementSets.phase_direction(ms)
    u,v,w = MeasurementSets.uvw(ms)
    ν = MeasurementSets.frequency(ms)
    ant1,ant2 = MeasurementSets.antennas(ms)

    frame = ReferenceFrame()
    set!(frame,MeasurementSets.position(ms))
    set!(frame,MeasurementSets.time(ms))
    sources = filter(source -> isabovehorizon(frame,source),sources)

    Nant  = numrows(Table(ms[kw"ANTENNA"]))
    Nfreq = length(ν)
    calibration = AmplitudeCalibration(Nant,Nfreq)

    data  = ms["DATA"]
    model = genvis(frame,phase_dir,sources,u,v,w,ν)
    flags = MeasurementSets.flags(ms)

    if force_imaging_columns || Tables.checkColumnExists(ms,"MODEL_DATA")
        ms["MODEL_DATA"] = model
    end

    ampcal!(calibration,data,model,flags,
            ant1,ant2,maxiter,tolerance)
    calibration
end

function ampcal!(calibration,data,model,flags,
                 ant1,ant2,maxiter,tolerance)
    for β = 1:Nfreq(calibration)
        # Calibrate the X polarization
        ampcal_onechannel!(slice(calibration.amplitudes,:,β,1),
                           slice(calibration.flags,:,β,1),
                           slice(data, 1,β,:),
                           slice(model,1,β,:),
                           slice(flags,1,β,:),
                           ant1, ant2,
                           maxiter, tolerance)
        # Calibrate the Y polarization
        ampcal_onechannel!(slice(calibration.amplitudes,:,β,2),
                           slice(calibration.flags,:,β,2),
                           slice(data, 4,β,:),
                           slice(model,4,β,:),
                           slice(flags,4,β,:),
                           ant1, ant2,
                           maxiter, tolerance)
    end
end

function ampcal_onechannel!(amplitudes, amplitude_flags,
                            data, model, data_flags,
                            ant1, ant2,
                            maxiter, tolerance)
    # If the entire channel is flagged, don't bother calibrating.
    if all(data_flags)
        amplitude_flags[:] = true
        return
    end

    square_data  = ampcal_makesquare( data,data_flags,ant1,ant2)
    square_model = ampcal_makesquare(model,data_flags,ant1,ant2)
    best_amplitudes = ones(Float64,length(amplitudes))
    converged = @iterate(AmpCalStep(),RK4,maxiter,tolerance,
                         best_amplitudes,square_data,square_model)

    # Flag the entire channel if the solution did not converge.
    if !converged
        amplitude_flags[:] = true
        return
    end

    # Output
    for ant = 1:length(amplitudes)
        amplitudes[ant] = best_amplitudes[ant]
    end
    nothing
end

function ampcal_makesquare(data,flags,ant1,ant2)
    gaincal_makesquare(data,flags,ant1,ant2)
end

"""
    ampcal_step(amplitudes,data,model) -> step

Given the `data` and `model` visibilities, and the current
guess for the gain `amplitudes`, solve for `step` such
that the new value of the amplitudes is `amplitudes+step`.
"""
function ampcal_step(input,data,model)
    Nant = length(input)
    step = zeros(Float64,Nant)
    @inbounds for j = 1:Nant
        numerator   = zero(Float64)
        denominator = zero(Float64)
        for i = 1:Nant
            GM = input[i]*model[i,j]
            numerator = numerator + abs(real(GM'*data[i,j]))
            denominator = denominator + abs2(GM)
        end
        ok = abs(denominator) > eps(Float64)
        step[j] = ifelse(ok,numerator/denominator - input[j],0)
    end
    step
end

immutable AmpCalStep <: StepFunction end
call(::AmpCalStep,g,V,M) = ampcal_step(g,V,M)

function corrupt!(data::Array{Complex64,3},
                  flags::Array{Bool,3},
                  cal::AmplitudeCalibration,
                  ant1,ant2)
    Nbase = length(ant1)
    for α = 1:Nbase, β = 1:Nfreq(cal)
        if (cal.flags[ant1[α],β,1] || cal.flags[ant1[α],β,2]
                                   || cal.flags[ant2[α],β,1]
                                   || cal.flags[ant2[α],β,2])
            flags[:,β,α] = true
        end
        data[1,β,α] = cal.amplitudes[ant1[α],β,1]*cal.amplitudes[ant2[α],β,1]*data[1,β,α]
        data[2,β,α] = cal.amplitudes[ant1[α],β,1]*cal.amplitudes[ant2[α],β,2]*data[2,β,α]
        data[3,β,α] = cal.amplitudes[ant1[α],β,2]*cal.amplitudes[ant2[α],β,1]*data[3,β,α]
        data[4,β,α] = cal.amplitudes[ant1[α],β,2]*cal.amplitudes[ant2[α],β,2]*data[4,β,α]
    end
end

