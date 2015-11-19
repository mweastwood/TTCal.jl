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

abstract Calibration

"""
    immutable GainCalibration <: Calibration

This type stores the information for calibrating the
electronic gains of the interferometer. That is, it stores
complex gains and flags for each antenna, frequency channel,
and polarization.
"""
immutable GainCalibration <: Calibration
    jones::Matrix{DiagonalJonesMatrix}
    flags::Matrix{Bool}
end

"""
    GainCalibration(Nant, Nfreq)

Create a calibration table for `Nant` antennas with
`Nfreq` frequency channels where all the gains are
initially set to unity.
"""
function GainCalibration(Nant,Nfreq)
    GainCalibration(ones(DiagonalJonesMatrix,Nant,Nfreq),zeros(Bool,Nant,Nfreq))
end

Base.similar(cal::GainCalibration) = GainCalibration(Nant(cal),Nfreq(cal))

"""
    immutable PolarizationCalibration <: Calibration

This type stores the information for calibrating the
polarization of the interferometer. That is, it stores
Jones matrices and flags for each antenna and each
frequency channel.
"""
immutable PolarizationCalibration <: Calibration
    jones::Matrix{JonesMatrix}
    flags::Matrix{Bool}
end

"""
    PolarizationCalibration(Nant, Nfreq)

Create a calibration table for `Nant` antennas with
`Nfreq` frequency channels where all the Jones matrices
are initially set to the identity matrix.
"""
function PolarizationCalibration(Nant,Nfreq)
    PolarizationCalibration(ones(JonesMatrix,Nant,Nfreq),zeros(Bool,Nant,Nfreq))
end

Base.similar(cal::PolarizationCalibration) = PolarizationCalibration(Nant(cal),Nfreq(cal))

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
    applycal!(ms::MeasurementSet, calibration::Calibration;
              apply_to_corrected = false, force_imaging_columns = false)

Apply the calibration to the given measurement set.

**Arguments:**

* `ms` - the measurement set to which the calibration will be applied
* `calibration` - the calibration that will be applied

**Keyword Arguments:**

* `apply_to_corrected` - if this is set to true, the calibration will be
    applied to the CORRECTED_DATA column instead of the DATA column
* `force_imaging_columns` - if this is set to true, the calibrated data
    will be written to the CORRECTED_DATA column regardless of whether
    or not the column already exists
"""
function applycal!(ms::MeasurementSet,
                   calibration::Calibration;
                   apply_to_corrected::Bool = false,
                   force_imaging_columns::Bool = false)
    data  = apply_to_corrected? get_corrected_data(ms) : get_data(ms)
    flags = get_flags(ms)
    applycal!(data,flags,calibration,ms.ant1,ms.ant2)
    set_corrected_data!(ms,data,force_imaging_columns)
    set_flags!(ms,flags)
    data
end

function applycal!(data::Array{Complex64,3},
                   flags::Array{Bool,3},
                   cal::Calibration,
                   ant1,ant2)
    inverse_cal = invert(cal)
    corrupt!(data,flags,inverse_cal,ant1,ant2)
    data
end

"""
    corrupt!(data::Array{Complex64,3}, cal::Calibration, ant1, ant2)

Corrupt the model data as if it had been observed
with an instrument with the given calibration.
"""
function corrupt!(data::Array{Complex64,3},
                  cal::Calibration,
                  ant1,ant2)
    flags = fill(false,size(data))
    corrupt!(data,flags,cal,ant1,ant2)
end

write(filename,calibration::Calibration) = JLD.save(File(format"JLD",filename),"cal",calibration)
read(filename) = JLD.load(filename,"cal")

function write_for_python(filename,calibration::GainCalibration)
    gains = zeros(Complex128,2,Nant(calibration),Nfreq(calibration))
    flags = zeros(      Bool,  Nant(calibration),Nfreq(calibration))
    for β = 1:Nfreq(calibration), ant = 1:Nant(calibration)
        gains[1,ant,β] = calibration.jones[ant,β].xx
        gains[2,ant,β] = calibration.jones[ant,β].yy
        flags[  ant,β] = calibration.flags[ant,β]
    end
    npzwrite(filename,Dict("gains" => gains, "flags" => flags))
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
    sources = abovehorizon(ms.frame,sources)
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
                maxiter, tolerance; quiet = false)
    quiet || (p = Progress(Nfreq(calibration), 1, "Thinking...", 50))
    for β = 1:Nfreq(calibration)
        solve_onechannel!(slice(calibration.jones,:,β),
                          slice(calibration.flags,:,β),
                          slice(data, :,β,:),
                          slice(model,:,β,:),
                          slice(flags,:,β,:),
                          ant1, ant2, maxiter, tolerance)
        quiet || next!(p)
    end
    if !quiet && sum(calibration.flags) > 0.5length(calibration.flags)
        warn("Over half of the calibration solutions are flagged.")
    end
end

function solve_onechannel!(jones, jones_flags,
                           data, model, data_flags,
                           ant1, ant2, maxiter, tolerance)
    # If the entire channel is flagged, don't bother calibrating.
    all(data_flags) && (jones_flags[:] = true; return)

    N = length(jones)
    T = eltype(jones)
    square_data  = makesquare( data,data_flags,ant1,ant2)
    square_model = makesquare(model,data_flags,ant1,ant2)
    best_jones   = copy(jones)

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
function stefcal_step{T}(input::AbstractVector{T},data,model)
    Nant = length(input)
    step = similar(input)
    @inbounds for j = 1:Nant
        numerator   = zero(T)
        denominator = zero(T)
        for i = 1:Nant
            GM = input[i]*model[i,j]
            V  = data[i,j]
            numerator   += inner_multiply(T,GM,V)
            denominator += inner_multiply(T,GM,GM)
        end
        ok = abs(det(denominator)) > eps(Float64)
        step[j] = ifelse(ok,(denominator\numerator)' - input[j],zero(T))
    end
    step
end

@inline function inner_multiply(::Type{DiagonalJonesMatrix},X,Y)
    DiagonalJonesMatrix(X.xx'*Y.xx + X.yx'*Y.yx,
                        X.xy'*Y.xy + X.yy'*Y.yy)
end

@inline inner_multiply(::Type{JonesMatrix},X,Y) = X'*Y

immutable StefCalStep <: StepFunction end
call(::StefCalStep,g,V,M) = stefcal_step(g,V,M)

