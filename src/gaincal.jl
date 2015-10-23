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

for (Cal,eltype) in ((:GainCalibration, Complex64),
                     (:AmplitudeCalibration, Float64))
    @eval immutable $Cal <: Calibration
        gains::Array{$eltype,3}
        flags::Array{Bool,3}
    end

    @eval function $Cal(Nant,Nfreq)
        gains =  ones($eltype,2,Nant,Nfreq)
        flags = zeros(   Bool,2,Nant,Nfreq)
        $Cal(gains,flags)
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
""" GainCalibration

@doc """
    immutable AmplitudeCalibration <: Calibration

This type stores the information for calibrating only
the amplitude of the electronic gains. It stores
the gain amplitudes and flags for each antenna, frequency channel,
and polarization.

    AmplitudeCalibration(Nant, Nfreq)

Create a calibration table for `Nant` antennas with
`Nfreq` frequency channels where all the amplitudes
are initially set to unity.
""" AmplitudeCalibration

"""
An alias for `Union{GainCalibration,AmplitudeCalibration}`.
This alis is useful because these calibration strategies
have very similar implementations.
"""
typealias ScalarCalibration Union{GainCalibration,AmplitudeCalibration}

Nant( cal::ScalarCalibration) = size(cal.gains,2)
Nfreq(cal::ScalarCalibration) = size(cal.gains,3)

"""
    corrupt!(data::Array{Complex64,3}, flags::Array{Bool,3},
             cal::ScalarCalibration, ant1, ant2)

Corrupt the data as if it was observed with the given calibration.
"""
function corrupt!(data::Array{Complex64,3}, flags::Array{Bool,3},
                  cal::ScalarCalibration, ant1, ant2)
    Nbase = length(ant1)
    for α = 1:Nbase, β = 1:Nfreq(cal)
        if (cal.flags[1,ant1[α],β] || cal.flags[2,ant1[α],β]
                                   || cal.flags[1,ant2[α],β]
                                   || cal.flags[2,ant2[α],β])
            flags[:,β,α] = true
        end
        data[1,β,α] = cal.gains[1,ant1[α],β]*conj(cal.gains[1,ant2[α],β])*data[1,β,α]
        data[2,β,α] = cal.gains[1,ant1[α],β]*conj(cal.gains[2,ant2[α],β])*data[2,β,α]
        data[3,β,α] = cal.gains[2,ant1[α],β]*conj(cal.gains[1,ant2[α],β])*data[3,β,α]
        data[4,β,α] = cal.gains[2,ant1[α],β]*conj(cal.gains[2,ant2[α],β])*data[4,β,α]
    end
end

doc"""
    invert(cal::ScalarCalibration)

Returns the inverse of the given calibration.
The gain $g$ of each antenna is set to $1/g$.
"""
function invert(cal::ScalarCalibration)
    output = similar(cal)
    for i in eachindex(output.gains,output.flags,
                          cal.gains,   cal.flags)
        output.gains[i] = inv(cal.gains[i])
        output.flags[i] = cal.flags[i]
    end
    output
end

"""
    fixphase!(cal::GainCalibration, reference_antenna)

Set the phase of the reference antenna and polarization to zero.

**Arguments:**

* `cal` - the calibration that will have its phase adjusted
* `reference_antenna` - a string containing the antenna number and polarization
    whose phase will be chosen to be zero (eg. "14y" or "62x")
"""
function fixphase!(cal::GainCalibration, reference_antenna)
    regex = r"(\d+)(x|y)"
    m = match(regex,reference_antenna)
    refant = parse(Int,m.captures[1])
    refpol = m.captures[2] == "x"? 1 : 2

    for β = 1:Nfreq(cal)
        factor = conj(cal.gains[refpol,refant,β]) / abs(cal.gains[refpol,refant,β])
        for ant = 1:Nant(cal), pol = 1:2
            cal.gains[pol,ant,β] = cal.gains[pol,ant,β]*factor
        end
    end
    cal
end

"""
    fixphase!(cal::AmplitudeCalibration, reference_antenna)

With an amplitude calibration there is no freedom to pick an arbitrary phase,
so this function does nothing.
"""
fixphase!(cal::AmplitudeCalibration, reference_antenna) = cal

"""
    gaincal(ms::MeasurementSet, sources::Vector{PointSource}, beam::BeamModel;
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
                 sources::Vector{PointSource},
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
           ms.ant1,ms.ant2,maxiter,tolerance,
           reference_antenna)
    calibration
end

"""
    ampcal(ms::MeasurementSet, sources::Vector{PointSource}, beam::BeamModel;
           maxiter = 20, tolerance = 1e-3, minuvw = 0.0,
           force_imaging_columns = false)

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
function ampcal(ms::MeasurementSet,
                sources::Vector{PointSource},
                beam::BeamModel;
                maxiter::Int = 20,
                tolerance::Float64 = 1e-3,
                minuvw::Float64 = 0.0,
                force_imaging_columns::Bool = false)
    sources = abovehorizon(ms.frame,sources)
    calibration = AmplitudeCalibration(ms.Nant,ms.Nfreq)
    data  = get_data(ms)
    model = genvis(ms,sources,beam)
    flags = get_flags(ms)
    flag_short_baselines!(flags,minuvw,ms.u,ms.v,ms.w,ms.ν)
    set_model_data!(ms,model,force_imaging_columns)
    solve!(calibration,data,model,flags,
           ms.ant1,ms.ant2,maxiter,tolerance)
    calibration
end

function solve!(calibration::ScalarCalibration,
                data, model, flags, ant1, ant2,
                maxiter, tolerance, reference_antenna)
    for β = 1:Nfreq(calibration)
        solve_scalar_onechannel!(slice(calibration.gains,:,:,β),
                                 slice(calibration.flags,:,:,β),
                                 slice(data, :,β,:),
                                 slice(model,:,β,:),
                                 slice(flags,:,β,:),
                                 ant1, ant2, maxiter, tolerance)
    end
    fixphase!(calibration,reference_antenna)
end

function solve_scalar_onechannel!(gains, gain_flags,
                                  data, model, data_flags,
                                  ant1, ant2, maxiter, tolerance)
    # If the entire channel is flagged, don't bother calibrating.
    all(data_flags) && (gain_flags[:] = true; return)

    N = length(gains)
    T = eltype(gains)
    square_data  = scalar_makesquare( data,data_flags,ant1,ant2)
    square_model = scalar_makesquare(model,data_flags,ant1,ant2)
    best_gains   = reshape(copy(gains),N)
    converged = @iterate(step(T),RK4,maxiter,tolerance,
                         best_gains,square_data,square_model)

    # Propagate antenna flags to the calibration solutions.
    # A flagged antenna should correspond to a row and column
    # of zeros in both square_data and square_model.
    antenna_flags = all(square_data .== 0,2)

    # Flag the entire channel if the solution did not converge.
    # However, we'll still write out our best guess for what
    # the gains should be.
    !converged && (antenna_flags[:] = true)

    gains[:] = best_gains
    gain_flags[:] = antenna_flags
    nothing
end

"""
    scalar_makesquare(data, flags, ant1, ant2)

Pack the data into a square Hermitian matrix such that
the data is ordered as follows:

    x₁x₁ x₁y₁ x₁x₂ x₁y₂
    y₁x₁ y₁y₁ y₁x₂ y₁y₂
    x₂x₁ y₂y₁ y₂x₂ y₂y₂
    y₂x₁ y₂y₁ y₂x₂ y₂y₂
                         .
                           .
                             .

Flagged correlations and autocorrelations are set to zero.
"""
function scalar_makesquare(data,flags,ant1,ant2)
    Nbase = length(ant1)
    Nant = maximum(ant1)
    output = zeros(Complex64,2Nant,2Nant)
    for α = 1:Nbase
        ant1[α] == ant2[α] && continue

        # indices into output for x/y polarizations and antennas 1/2
        x1 = 2ant1[α] - 1
        y1 = 2ant1[α] - 0
        x2 = 2ant2[α] - 1
        y2 = 2ant2[α] - 0

        output[x1,x2] = ifelse(flags[1,α],0,data[1,α])
        output[x1,y2] = ifelse(flags[2,α],0,data[2,α])
        output[y1,x2] = ifelse(flags[3,α],0,data[3,α])
        output[y1,y2] = ifelse(flags[4,α],0,data[4,α])

        output[x2,x1] = conj(output[x1,x2])
        output[y2,x1] = conj(output[x1,y2])
        output[x2,y1] = conj(output[y1,x2])
        output[y2,y1] = conj(output[y1,y2])
    end
    output
end

doc"""
    gaincal_step(gains,data,model) -> step

Given the `data` and `model` visibilities, and the current
guess for the electronic `gains`, solve for `step` such
that the new value of the gains is `gains+step`.

The update step is defined such that the new value of the
gains minimizes

\[
    \sum_{i,j}\|V_{i,j} - g_i g_{j,new}^* M_{i,j}\|^2$,
\]

where $i$ and $j$ label the antennas, $V$ labels the measured
visibilities, $M$ labels the model visibilities, and $g$
labels the complex gains.

*References:*

* Michell, D. et al. 2008, JSTSP, 2, 5.
* Salvini, S. & Wijnholds, S. 2014, A&A, 571, 97.
"""
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
        end
        ok = abs(denominator) > eps(Float32)
        step[j] = ifelse(ok,conj(numerator/denominator) - input[j],0)
    end
    step
end

immutable GainCalStep <: StepFunction end
call(::GainCalStep,g,V,M) = gaincal_step(g,V,M)
step(::Type{Complex64})   = GainCalStep()

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
step(::Type{Float64})    = AmpCalStep()

