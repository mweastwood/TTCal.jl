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
    immutable PolarizationCalibration <: Calibration

This type stores the information for calibrating the
polarization of the interferometer. That is, it stores
Jones matrices and flags for each antenna and each
frequency channel.
"""
immutable PolarizationCalibration <: Calibration
    jones::Array{JonesMatrix,2}
    flags::Array{Bool,2}
end

"""
    PolarizationCalibration(Nant, Nfreq)

Create a calibration table for `Nant` antennas with
`Nfreq` frequency channels where all the Jones matrices
are initially set to the identity matrix.
"""
function PolarizationCalibration(Nant,Nfreq)
    jones = fill(JonesMatrix(),Nant,Nfreq)
    flags = fill(false,Nant,Nfreq)
    PolarizationCalibration(jones,flags)
end

Nant( cal::PolarizationCalibration) = size(cal.jones,1)
Nfreq(cal::PolarizationCalibration) = size(cal.jones,2)

"""
    corrupt!(data::Array{Complex64,3}, flags::Array{Bool,3},
             cal::PolarizationCalibration, ant1, ant2)

Corrupt the data as if it was observed with the given calibration.
"""
function corrupt!(data::Array{Complex64,3}, flags::Array{Bool,3},
                  cal::PolarizationCalibration, ant1, ant2)
    Nbase = length(ant1)
    for α = 1:Nbase, β = 1:Nfreq(cal)
        if cal.flags[ant1[α],β] || cal.flags[ant2[α],β]
            flags[:,β,α] = true
        end
        V = JonesMatrix(data[1,β,α],data[2,β,α],
                        data[3,β,α],data[4,β,α])
        G1 = cal.jones[ant1[α],β]
        G2 = cal.jones[ant2[α],β]
        V = G1*V*G2'
        data[1,β,α] = V.xx
        data[2,β,α] = V.xy
        data[3,β,α] = V.yx
        data[4,β,α] = V.yy
    end
end

doc"""
    invert(cal::PolarizationCalibration)

Returns the inverse of the given calibration.
The Jones matrix $J$ of each antenna is set to $J^{-1}$.
"""
function invert(cal::PolarizationCalibration)
    output = PolarizationCalibration(Nant(cal),Nfreq(cal))
    for i in eachindex(output.jones,output.flags,cal.jones,cal.flags)
        output.jones[i] = inv(cal.jones[i])
        output.flags[i] = cal.flags[i]
    end
    output
end

"""
    polcal(ms::MeasurementSet, sources::Vector{PointSources}, beam::BeamModel;
           maxiter = 20, tolerance = 1e-3, minuvw = 0.0,
           force_imaging_columns = false)

Solve for the polarization properties of the interferometer.

**Arguments:**

* `ms` - the measurement set from which to derive the calibration
* `sources` - the list of points sources to use as the sky model
* `beam` - the beam model

**Keyword Arguments:**

* `maxiter` - the maximum number of Runge-Kutta steps to take on each
    frequency channel
* `tolerance` - the relative tolerance to use while checking to see if
    more iterations are required
* `minuvw` - the minimum baseline length (measured in wavelengths) to be
    used during the calibration procedure
* `force_imaging_columns` - if this is set to true, the MODEL_DATA column
    will be created and populated with model visibilities even if it
    doesn't already exist
"""
function polcal(ms::MeasurementSet,
                sources::Vector{PointSource},
                beam::BeamModel;
                maxiter::Int = 20,
                tolerance::Float64 = 1e-3,
                minuvw::Float64 = 0.0,
                force_imaging_columns::Bool = false)
    sources = filter(source -> isabovehorizon(ms.frame,source),sources)
    calibration = PolarizationCalibration(ms.Nant,ms.Nfreq)
    data  = get_corrected_data(ms)
    model = genvis(ms,sources,beam)
    flags = get_flags(ms)
    flag_short_baselines!(flags,minuvw,ms.u,ms.v,ms.w,ms.ν)
    set_model_data!(ms,model)
    polcal!(calibration,data,model,flags,
            ms.ant1,ms.ant2,maxiter,tolerance)
    calibration
end

function solve!(calibration::PolarizationCalibration,
                data, model, flags, ant1, ant2,
                maxiter, tolerance)
    for β = 1:Nfreq(calibration)
        polcal_onechannel!(slice(calibration.jones,:,β),
                           slice(calibration.flags,:,β),
                           slice( data,:,β,:),
                           slice(model,:,β,:),
                           slice(flags,:,β,:),
                           ant1, ant2,
                           maxiter, tolerance)
    end
    # TODO: In the unpolarized gain calibration, we are free to
    # pick the overall phase. Usually this is used to zero the
    # phase of a reference antenna. Is there an equivalent
    # degree of freedom for the polarized case?
end

function solve_jones_onechannel!(jones, jones_flags,
                                 data, model, data_flags,
                                 ant1, ant2, maxiter, tolerance)
    # If the entire channel is flagged, don't bother calibrating.
    all(data_flags) && (gain_flags[:] = true; return)

    square_data  = polcal_makesquare( data,data_flags,ant1,ant2)
    square_model = polcal_makesquare(model,data_flags,ant1,ant2)
    best_jones   = copy(jones)
    converged = @iterate(PolCalStep(),RK4,maxiter,tolerance,
                         best_jones,square_data,square_model)

    # Propagate antenna flags to the calibration solutions.
    # A flagged antenna should correspond to a row and column
    # of zeros in both square_data and square_model.
    antenna_flags = all(square_data .== zero(JonesMatrix),2)

    # Flag the entire channel if the solution did not converge.
    # However, we'll still write out our best guess for what
    # the jones matrices should be.
    !converged && (antenna_flags[:] = true)

    jones[:] = best_jones
    jones_flags[:] = antenna_flags
    nothing
end

"""
    polcal_makesquare(data, flags, ant1, ant2)

Pack the data into a square Hermitian matrix where each element
is a Jones matrix.

Compare this to the regular complex gain calibration where each
element is a complex scalar.

The packing order is:

    11 12 13
    21 22 23
    31 32 33
             .
               .
                 .
"""
function polcal_makesquare(data,flags,ant1,ant2)
    Nbase = size(data,2)
    Nant  = round(Integer,div(sqrt(1+8Nbase)-1,2))
    output = zeros(JonesMatrix,Nant,Nant)
    for α = 1:Nbase
        (flags[1,α] || flags[2,α] || flags[3,α] || flags[4,α]) && continue
        ant1[α] == ant2[α] && continue
        output[ant1[α],ant2[α]] = JonesMatrix(data[1,α],data[2,α],
                                              data[3,α],data[4,α])
        output[ant2[α],ant1[α]] = output[ant1[α],ant2[α]]'
    end
    output
end

"""
    polcal_step(jones,data,model) -> step

Given the `data` and `model` visibilities, and the current
guess for the Jones matrices (`jones`), solve for `step` such
that the new value of the Jones matrices is `jones+step`.
"""
function polcal_step(input,data,model)
    Nant = length(input)
    step = fill(JonesMatrix(),Nant)
    @inbounds for j = 1:Nant
        numerator   = zero(JonesMatrix)
        denominator = zero(JonesMatrix)
        for i = 1:Nant
            GM = input[i]*model[i,j]
            numerator = numerator + GM'*data[i,j]
            denominator = denominator + GM'*GM
        end
        ok = abs(det(denominator)) > eps(Float32)
        step[j] = ifelse(ok,conj(denominator\numerator) - input[j],zero(JonesMatrix))
    end
    step
end

immutable PolCalStep <: StepFunction end
call(::PolCalStep,input,data,model) = polcal_step(input,data,model)

