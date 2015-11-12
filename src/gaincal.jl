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

for (Cal,eltype) in ((:GainCalibration, DiagonalJonesMatrix),
                     (:PolarizationCalibration, JonesMatrix))
    @eval immutable $Cal <: Calibration
        jones::Array{$eltype,2}
        flags::Array{Bool,2}
    end

    @eval function $Cal(Nant,Nfreq)
        jones =  ones($eltype,Nant,Nfreq)
        flags = zeros(   Bool,Nant,Nfreq)
        $Cal(jones,flags)
    end

    @eval Base.similar(cal::$Cal) = $Cal(Nant(cal),Nfreq(cal))
end

@doc """
    immutable GainCalibration <: Calibration

This type stores the information for calibrating the
electronic gains of the interferometer. That is, it stores
complex gains and flags for each antenna, frequency channel,
and polarization.

    GainCalibration(Nant, Nfreq)

Create a calibration table for `Nant` antennas with
`Nfreq` frequency channels where all the gains are
initially set to unity.
""" ->
GainCalibration

@doc """
    immutable PolarizationCalibration <: Calibration

This type stores the information for calibrating the
polarization of the interferometer. That is, it stores
Jones matrices and flags for each antenna and each
frequency channel.

    PolarizationCalibration(Nant, Nfreq)

Create a calibration table for `Nant` antennas with
`Nfreq` frequency channels where all the Jones matrices
are initially set to the identity matrix.
""" ->
PolarizationCalibration

Nant( cal::Calibration) = size(cal.jones,1)
Nfreq(cal::Calibration) = size(cal.jones,2)

"""
    corrupt!(data::Array{Complex64,3}, flags::Array{Bool,3},
             cal::Calibration, ant1, ant2)

Corrupt the data as if it was observed with the given calibration.
"""
function corrupt!(data::Array{Complex64,3}, flags::Array{Bool,3},
                  cal::Calibration, ant1, ant2)
    Nbase = length(ant1)
    for α = 1:Nbase, β = 1:Nfreq(cal)
        if cal.flags[ant1[α],β] || cal.flags[ant2[α],β]
            flags[:,β,α] = true
        end
        V = JonesMatrix(data[1,β,α],data[2,β,α],
                        data[3,β,α],data[4,β,α])
        J1 = cal.jones[ant1[α],β]
        J2 = cal.jones[ant2[α],β]
        V = J1*V*J2'
        data[1,β,α] = V.xx
        data[2,β,α] = V.xy
        data[3,β,α] = V.yx
        data[4,β,α] = V.yy
    end
end

doc"""
    invert(cal::Calibration)

Returns the inverse of the given calibration.
The Jones matrix $J$ of each antenna is set to $J^{-1}$.
"""
function invert(cal::Calibration)
    output = similar(cal)
    for i in eachindex(output.jones,output.flags,
                          cal.jones,   cal.flags)
        output.jones[i] = inv(cal.jones[i])
        output.flags[i] = cal.flags[i]
    end
    output
end

"""
    fixphase!(cal::Calibration, reference_antenna)

Set the phase of the reference antenna and polarization to zero.

**Arguments:**

* `cal` - the calibration that will have its phase adjusted
* `reference_antenna` - a string containing the antenna number and polarization
    whose phase will be chosen to be zero (eg. "14y" or "62x")
"""
function fixphase!(cal::Calibration, reference_antenna)
    regex = r"(\d+)(x|y)"
    m = match(regex,reference_antenna)
    refant = parse(Int,m.captures[1])
    refpol = m.captures[2] == "x"? 1 : 2

    for β = 1:Nfreq(cal)
        ref = refpol == 1? cal.jones[refant,β].xx : cal.jones[refant,β].yy
        factor = conj(ref) / abs(ref)
        for ant = 1:Nant(cal)
            cal.jones[ant,β] = cal.jones[ant,β]*factor
        end
    end
    cal
end

"""
    gaincal(ms::MeasurementSet, sources::Vector{Source}, beam::BeamModel;
            maxiter = 20, tolerance = 1e-3, minuvw = 0.0,
            reference_antenna = "1x", force_imaging_columns = false)

Solve for the interferometer's electronic gains.

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
* `reference_antenna` - a string containing the antenna number and polarization
    whose phase will be chosen to be zero (eg. "14y" or "62x")
* `force_imaging_columns` - if this is set to true, the MODEL_DATA column
    will be created and populated with model visibilities even if it
    doesn't already exist
"""
function gaincal(ms::MeasurementSet,
                 sources::Vector{Source},
                 beam::BeamModel;
                 maxiter::Int = 20,
                 tolerance::Float64 = 1e-3,
                 minuvw::Float64 = 0.0,
                 reference_antenna::ASCIIString = "1x",
                 force_imaging_columns::Bool = false)
    sources = abovehorizon(ms.frame,sources)
    calibration = GainCalibration(ms.Nant,ms.Nfreq)
    data  = get_data(ms)
    model = genvis(ms,sources,beam)
    flags = get_flags(ms)
    flag_short_baselines!(flags,minuvw,ms.u,ms.v,ms.w,ms.ν)
    set_model_data!(ms,model,force_imaging_columns)
    solve!(calibration,data,model,flags,
           ms.ant1,ms.ant2,maxiter,tolerance)
    fixphase!(calibration,reference_antenna)
    calibration
end

"""
    polcal(ms::MeasurementSet, sources::Vector{Source}, beam::BeamModel;
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
                sources::Vector{Source},
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

function solve!(calibration::Calibration,
                data, model, flags, ant1, ant2,
                maxiter, tolerance)
    for β = 1:Nfreq(calibration)
        solve_onechannel!(slice(calibration.jones,:,β),
                          slice(calibration.flags,:,β),
                          slice(data, :,β,:),
                          slice(model,:,β,:),
                          slice(flags,:,β,:),
                          ant1, ant2, maxiter, tolerance)
    end
end

function solve_onechannel!(jones, jones_flags,
                           data, model, data_flags,
                           ant1, ant2, maxiter, tolerance)
    # If the entire channel is flagged, don't bother calibrating.
    all(data_flags) && (gain_flags[:] = true; return)

    N = length(jones)
    T = eltype(jones)
    square_data  = makesquare( data,data_flags,ant1,ant2)
    square_model = makesquare(model,data_flags,ant1,ant2)
    best_jones   = copy(jones)

    #step = gaincal_step(ones(Complex64,256),gaincal_makesquare(data[1,:],data_flags[1,:],ant1,ant2),
    #                                        gaincal_makesquare(model[1,:],data_flags[1,:],ant1,ant2))
    #@show step[1]
    #println("-------")
#
#    step = gaincal_step(ones(Complex64,256),gaincal_makesquare(data[4,:],data_flags[4,:],ant1,ant2),
#                                            gaincal_makesquare(model[4,:],data_flags[4,:],ant1,ant2))
#    @show step[1]
#    println("-------")
#
#    step = stefcal_step(best_jones,square_data,square_model)
#    @show step[1]
#    #@time stefcal_step(best_jones,square_data,square_model)
#    #println(@code_llvm stefcal_step(best_jones,square_data,square_model))
#    error("stop")

    converged = @iterate(StefCalStep(),RK4,maxiter,tolerance,
                         best_jones,square_data,square_model)

    # Propagate antenna flags to the calibration solutions.
    # A flagged antenna should correspond to a row and column
    # of zeros in both square_data and square_model.
    antenna_flags = all(square_data .== zero(JonesMatrix),2)

    # Flag the entire channel if the solution did not converge.
    # However, we'll still write out our best guess for what
    # the gains should be.
    converged || (antenna_flags[:] = true)

    jones[:] = best_jones
    jones_flags[:] = antenna_flags
    jones,jones_flags
end

doc"""
    makesquare(data, flags, ant1, ant2)

Pack the data into a square Hermitian matrix such that
the data is ordered as follows:

\\[
    \begin{pmatrix}
        V_{11} & V_{12} & V_{13} &        & \\\\
        V_{21} & V_{22} & V_{23} &        & \\\\
        V_{31} & V_{32} & V_{33} &        & \\\\
               &        &        & \ddots & \\\\
    \end{pmatrix}
\\]

Flagged correlations and autocorrelations are set to zero.
"""
function makesquare(data, flags, ant1, ant2)
    Nbase = length(ant1)
    Nant = maximum(ant1)
    output = zeros(JonesMatrix,Nant,Nant)
    for α = 1:Nbase
        (flags[1,α] || flags[2,α] || flags[3,α] || flags[4,α]) && continue
        ant1[α] == ant2[α] && continue

        V = JonesMatrix(data[1,α],data[2,α],
                        data[3,α],data[4,α])
        #V = JonesMatrix(data[1,α],0,
        #                0,data[4,α])
        output[ant1[α],ant2[α]] = V
        output[ant2[α],ant1[α]] = output[ant1[α],ant2[α]]'
    end
    output
end

doc"""
    stefcal_step(input,data,model) -> step

Given the `data` and `model` visibilities, and the current
guess for the Jones matrices, solve for `step` such
that the new value of the Jones matrices is `input+step`.

The update step is defined such that the new value of the
Jones matrices minimizes

\\[
    \sum_{i,j}\|V_{i,j} - G_i M_{i,j} G_{j,new}^*\|^2$,
\\]

where $i$ and $j$ label the antennas, $V$ labels the measured
visibilities, $M$ labels the model visibilities, and $G$
labels the Jones matrices.

*References:*

* Michell, D. et al. 2008, JSTSP, 2, 5.
* Salvini, S. & Wijnholds, S. 2014, A&A, 571, 97.
"""
function stefcal_step(input,data,model)
    T = eltype(input)
    Nant = length(input)
    step = similar(input)
    @inbounds for j = 1:Nant
        numerator   = zero(T)
        denominator = zero(T)
        for i = 1:Nant
            GM = input[i]*model[i,j]
            V  = data[i,j]
            # note that the compiler should remove the dead branch
            #if T == DiagonalJonesMatrix
                numerator   += DiagonalJonesMatrix(GM.xx'*V.xx + GM.yx'*V.yx,
                                                   GM.xy'*V.xy + GM.yy'*V.yy)
                denominator += DiagonalJonesMatrix(GM.xx'*GM.xx + GM.yx'*GM.yx,
                                                   GM.xy'*GM.xy + GM.yy'*GM.yy)
            #else
            #    numerator   += GM'*V
            #    denominator += GM'*GM
            #end
        end
        ok = abs(det(denominator)) > eps(Float64)
        step[j] = ifelse(ok,(denominator\numerator)' - input[j],zero(T))
    end
    step
end

immutable StefCalStep <: StepFunction end
call(::StefCalStep,g,V,M) = stefcal_step(g,V,M)

#=
function gaincal_step(input,data,model)
    Nant = length(input)
    step = zeros(Complex64,Nant)
    @inbounds for j = 1:Nant
        numerator   = zero(Complex64)
        denominator = zero(Complex64)
        for i = 1:Nant
            GM = input[i]*model[i,j]
            numerator = numerator + GM'*data[i,j]
            denominator = denominator + GM'*GM
            if i == 256 && j == 1
                @show GM GM'*data[i,j] GM'*GM
                @show numerator
                @show denominator
            end
        end
        ok = abs(denominator) > eps(Float32)
        step[j] = ifelse(ok,conj(numerator/denominator) - input[j],0)
    end
    step
end

function gaincal_makesquare(data,flags,ant1,ant2)
    Nbase = length(data)
    Nant  = round(Integer,div(sqrt(1+8Nbase)-1,2))
    output = zeros(Complex64,Nant,Nant)
    for α = 1:Nbase
        flags[α] && continue
        ant1[α] == ant2[α] && continue
        output[ant1[α],ant2[α]] = data[α]
        output[ant2[α],ant1[α]] = conj(data[α])
    end
    output
end
=#

