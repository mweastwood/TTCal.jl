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

macro generate_calibration(name, jones)
    quote
        Base.@__doc__ immutable $name <: Calibration
            jones :: Matrix{$jones}
            flags :: Matrix{Bool}
        end
        function $name(Nant, Nfreq)
            $name(ones($jones, Nant, Nfreq), zeros(Bool, Nant, Nfreq))
        end
        Base.similar(cal::$name) = $name(Nant(cal), Nfreq(cal))
    end |> esc
end

"""
    GainCalibration

*Description*

This type represents the gain calibration of an interferometer.

Each antenna and frequency channel receives a diagonal Jones matrix.
The diagonal terms of the Jones matrix are the complex gains of the `x`
and `y` polarization respectively. The off-diagonal terms represent the
polarization leakage from `x` to `y` and `y` to `x`. These off-diagonal
terms are assumed to be zero.

The `gaincal` routine is used to solve for the interferometer's
gain calibration.

*Fields*

* `jones` - an array of diagonal Jones matrices (one per antenna and frequency channel)
* `flags` - a corresponding list of flags
"""
@generate_calibration GainCalibration DiagonalJonesMatrix

"""
    PolarizationCalibration

*Description*

This type represents the polarization calibration of an interferometer.

Each antenna and frequency channel receives a full Jones matrix.
The diagonal terms of the Jones matrix are the complex gains of the `x`
and `y` polarization respectively. The off-diagonal terms represent the
polarization leakage from `x` to `y` and `y` to `x`. All of these
terms are included in this calibration.

The `polcal` routine is used to solve for the interferometer's
polarization calibration.

*Fields*

* `jones` - an array of Jones matrices (one per antenna and frequency channel)
* `flags` - a corresponding list of flags
"""
@generate_calibration PolarizationCalibration JonesMatrix

Nant( cal::Calibration) = size(cal.jones, 1)
Nfreq(cal::Calibration) = size(cal.jones, 2)

write(filename, calibration::Calibration) = JLD.save(File(format"JLD", filename), "cal", calibration)
read(filename) = JLD.load(filename, "cal")

function write_for_python(filename, calibration::GainCalibration)
    gains = zeros(Complex128, 2, Nant(calibration), Nfreq(calibration))
    flags = zeros(      Bool,    Nant(calibration), Nfreq(calibration))
    for β = 1:Nfreq(calibration), ant = 1:Nant(calibration)
        gains[1,ant,β] = calibration.jones[ant,β].xx
        gains[2,ant,β] = calibration.jones[ant,β].yy
        flags[  ant,β] = calibration.flags[ant,β]
    end
    npzwrite(filename, Dict("gains" => gains, "flags" => flags))
end

function write_for_python(filename, calibration::PolarizationCalibration)
    gains = zeros(Complex128, 4, Nant(calibration), Nfreq(calibration))
    flags = zeros(      Bool,    Nant(calibration), Nfreq(calibration))
    for β = 1:Nfreq(calibration), ant = 1:Nant(calibration)
        gains[1,ant,β] = calibration.jones[ant,β].xx
        gains[2,ant,β] = calibration.jones[ant,β].xy
        gains[3,ant,β] = calibration.jones[ant,β].yx
        gains[4,ant,β] = calibration.jones[ant,β].yy
        flags[  ant,β] = calibration.flags[ant,β]
    end
    npzwrite(filename, Dict("gains" => gains, "flags" => flags))
end

# corrupt / applycal

doc"""
    corrupt!(visibilities, metadata, calibration)

*Description*

Corrupt the visibilities as if they were observed with the given calibration.
That is if we have the model visibility $V_{i,j}$ on baseline $i,j$, and the
corresponding Jones matrices $J_i$ and $J_j$ for antennas $i$ and $j$ compute

$$V_{i,j} \rightarrow J_i V_{i,j} J_j^*$$

*Arguments*

* `visibilities` - the list of visibilities to corrupt
* `metadata` - the metadata describing the interferometer
* `calibration` - the calibration in question
"""
function corrupt!(visibilities::Visibilities, meta::Metadata, calibration::Calibration)
    for β = 1:Nfreq(meta), α = 1:Nbase(meta)
        antenna1 = meta.baselines[α].antenna1
        antenna2 = meta.baselines[α].antenna2
        if calibration.flags[antenna1,β] || calibration.flags[antenna2,β]
            visibilities.flags[α,β] = true
        end
        V  = visibilities.data[α,β]
        J₁ = calibration.jones[antenna1,β]
        J₂ = calibration.jones[antenna2,β]
        visibilities.data[α,β] = J₁*V*J₂'
    end
    visibilities
end


doc"""
    applycal!(visibilities, metadata, calibration)

*Description*

Apply the calibration to the given visibilities.
That is if we have the model visibility $V_{i,j}$ on baseline $i,j$, and the
corresponding Jones matrices $J_i$ and $J_j$ for antennas $i$ and $j$ compute

$$V_{i,j} \rightarrow J_i^{-1} V_{i,j} (J_j^{-1})^*$$

*Arguments*

* `visibilities` - the list of visibilities to corrupt
* `metadata` - the metadata describing the interferometer
* `calibration` - the calibration in question
"""
function applycal!(visibilities::Visibilities, meta::Metadata, calibration::Calibration)
    inverse_cal = invert(calibration)
    corrupt!(visibilities, meta, inverse_cal)
    visibilities
end

doc"""
    invert(calibration)

Returns the inverse of the given calibration.
The Jones matrices $J$ of each antenna is set to $J^{-1}$.
"""
function invert(calibration::Calibration)
    output = similar(calibration)
    for i in eachindex(     output.jones,      output.flags,
                       calibration.jones, calibration.flags)
        output.jones[i] = inv(calibration.jones[i])
        output.flags[i] = calibration.flags[i]
    end
    output
end

# gaincal / polcal

"""
    gaincal(visibilities, metadata, sources)
    gaincal(visibilities, metadata, model_visibilities)

*Description*

Solve for the interferometer's electronic gains.

*Arguments*

* `visibilities` - the visibilities measured by the interferometer
* `metadata` - the metadata describing the interferometer
* `sources` - a list of sources comprising the sky model
* `model_visibilities` - alternatively the sky model visibilities can be provided

*Keyword Arguments*

* `maxiter` - the maximum number of iterations to take on each frequency channel (defaults to `20`)
* `tolerance` - the relative tolerance used to test for convergence (defaults to `1e-3`)
* `quiet` - suppresses printing if set to `true` (defaults to `false`)
"""
function gaincal{S<:Source}(visibilities::Visibilities, meta::Metadata, sources::Vector{S};
                            maxiter = 20, tolerance = 1e-3, quiet = false)
    frame = reference_frame(meta)
    sources = abovehorizon(frame, sources)
    model = genvis(meta, sources)
    gaincal(visibilities, meta, model, maxiter=maxiter, tolerance=tolerance, quiet=quiet)
end

function gaincal(visibilities::Visibilities, meta::Metadata, model::Visibilities;
                 maxiter = 20, tolerance = 1e-3, quiet = false)
    calibration = GainCalibration(Nant(meta), Nfreq(meta))
    solve!(calibration, visibilities, model, meta, maxiter, tolerance, quiet)
    calibration
end

"""
    polcal(visibilities, metadata, sources)
    polcal(visibilities, metadata, model_visibilities)

*Description*

Solve for the polarization properties of the interferometer.

*Arguments*

* `visibilities` - the visibilities measured by the interferometer
* `metadata` - the metadata describing the interferometer
* `sources` - a list of sources comprising the sky model
* `model_visibilities` - alternatively the sky model visibilities can be provided

*Keyword Arguments*

* `maxiter` - the maximum number of iterations to take on each frequency channel (defaults to `20`)
* `tolerance` - the relative tolerance used to test for convergence (defaults to `1e-3`)
* `quiet` - suppresses printing if set to `true` (defaults to `false`)
"""
function polcal{S<:Source}(visibilities::Visibilities, meta::Metadata, sources::Vector{S};
                           maxiter::Int = 20, tolerance::Float64 = 1e-3, quiet = false)
    frame = reference_frame(meta)
    sources = abovehorizon(frame, sources)
    model = genvis(meta, sources)
    polcal(visibilities, meta, model, maxiter=maxiter, tolerance=tolerance, quiet=quiet)
end

function polcal(visibilities::Visibilities, meta::Metadata, model::Visibilities;
                maxiter = 20, tolerance = 1e-3, quiet = false)
    calibration = PolarizationCalibration(Nant(meta), Nfreq(meta))
    solve!(calibration, visibilities, model, meta, maxiter, tolerance, quiet)
    calibration
end

function solve!(calibration, measured_visibilities, model_visibilities, metadata, maxiter, tolerance, quiet)
    square_measured, square_model = makesquare(measured_visibilities, model_visibilities, metadata)
    quiet || (p = Progress(Nfreq(metadata), "Calibrating: "))
    for β = 1:Nfreq(metadata)
        solve_onechannel!(slice(calibration.jones, :, β),
                          slice(calibration.flags, :, β),
                          slice(square_measured, :, :, β),
                          slice(square_model,    :, :, β),
                          maxiter, tolerance)
        quiet || next!(p)
    end
    if !quiet
        # Print a summary of the flags
        flagged_channels = sum(all(calibration.flags, 1))
        flagged_antennas = sum(all(calibration.flags, 2))
        percentage = 100 * sum(calibration.flags) / length(calibration.flags)
        @printf("(%d antennas flagged, %d channels flagged, %0.2f percent total)\n",
                flagged_antennas, flagged_channels, percentage)
    end
end

"Solve one channel at a time."
function solve_onechannel!(jones, flags, measured, model, maxiter, tolerance)
    converged = iterate(stefcalstep, RK4, maxiter, tolerance, false, jones, measured, model)
    flag_solution!(jones, flags, measured, model, converged)
    jones
end

"Solve with one solution for all channels."
function solve_allchannels!(calibration, measured, model, metadata, maxiter, tolerance)
    G = slice(calibration.jones, :, 1)
    F = fill(false, length(G))
    V, M = makesquare(measured, model, metadata)
    converged = iterate(stefcalstep, RK4, maxiter, tolerance, false, G, V, M)
    flag_solution!(G, F, V, M, converged)
    for β = 1:size(calibration.jones, 2)
        calibration.jones[:,β] = G
        calibration.flags[:,β] = F
    end
end

doc"""
Pack the visibilities into a square Hermitian matrices such that
the visibilities are ordered as follows:

$$\begin{pmatrix}
V_{11} & V_{12} & V_{13} &        & \\\\
V_{21} & V_{22} & V_{23} &        & \\\\
V_{31} & V_{32} & V_{33} &        & \\\\
       &        &        & \ddots & \\\\
\end{pmatrix}$$

Auto-correlations and flagged cross-correlations are set to zero.
"""
function makesquare(measured, model, meta)
    output_measured = zeros(JonesMatrix, Nant(meta), Nant(meta), Nfreq(meta))
    output_model    = zeros(JonesMatrix, Nant(meta), Nant(meta), Nfreq(meta))
    for β = 1:Nfreq(meta), α = 1:Nbase(meta)
        measured.flags[α,β] && continue
        antenna1 = meta.baselines[α].antenna1
        antenna2 = meta.baselines[α].antenna2
        antenna1 == antenna2 && continue
        output_measured[antenna1,antenna2,β] = measured.data[α,β]
        output_measured[antenna2,antenna1,β] = measured.data[α,β]'
        output_model[antenna1,antenna2,β] = model.data[α,β]
        output_model[antenna2,antenna1,β] = model.data[α,β]'
    end
    output_measured, output_model
end

"""
Inspect the calibration solution and apply flags as necessary.
"""
function flag_solution!(jones, flags, measured, model, converged)
    Nant  = length(jones)

    # Flag everything if the solution did not converge.
    if !converged
        flags[:] = true
        return flags
    end

    # Find flagged antennas by looking for columns of zeros.
    for j = 1:Nant
        isflagged = true
        for i = 1:Nant
            if measured[i,j] != zero(JonesMatrix)
                isflagged = false
                break
            end
        end
        if isflagged
            flags[j] = true
        end
    end

    flags
end

"""
    fixphase!(calibration, reference_antenna)

**Description**

Set the phase of the reference antenna and polarization to zero.

*Arguments*

* `calibration` - the calibration that will have its phase adjusted
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

# step functions

doc"""
    stefcal_step(input, measured, model)

**Description**

Given the `measured` and `model` visibilities, and the current
guess for the Jones matrices, solve for `step` such
that the new value of the Jones matrices is `input+step`.

The update step is defined such that the new value of the
Jones matrices minimizes

$$\sum_{i,j}\|V_{i,j} - J_i M_{i,j} J_{j,\rm new}^*\|^2,$$

where $i$ and $j$ label the antennas, $V$ labels the measured
visibilities, $M$ labels the model visibilities, and $J$
labels the Jones matrices.

*References*

* Michell, D. et al. 2008, JSTSP, 2, 5.
* Salvini, S. & Wijnholds, S. 2014, A&A, 571, 97.
"""
function stefcal_step{T}(input::AbstractVector{T}, measured::Matrix, model::Matrix)
    Nant = length(input) # number of antennas
    step = similar(input)
    @inbounds for j = 1:Nant
        numerator   = zero(T)
        denominator = zero(T)
        for i = 1:Nant
            GM = input[i]*model[i,j]
            V  = measured[i,j]
            numerator   += inner_multiply(T, GM, V)
            denominator += inner_multiply(T, GM, GM)
        end
        ok = abs(det(denominator)) > eps(Float64)
        step[j] = ifelse(ok, (denominator\numerator)' - input[j], zero(T))
    end
    step
end

function stefcal_step{T}(input::AbstractVector{T}, measured, model)
    # This version is called if `measured` and `model` have 3 dimensions
    # indicating that we want to use multiple integrations or multiple
    # frequency channels to solve for the calibration.
    Nant = length(input)     # number of antennas
    Nint = size(measured, 3) # number of integrations
    step = similar(input)
    @inbounds for j = 1:Nant
        numerator   = zero(T)
        denominator = zero(T)
        for t = 1:Nint, i = 1:Nant
            GM = input[i]*model[i,j,t]
            V  = measured[i,j,t]
            numerator   += inner_multiply(T, GM, V)
            denominator += inner_multiply(T, GM, GM)
        end
        ok = abs(det(denominator)) > eps(Float64)
        step[j] = ifelse(ok, (denominator\numerator)' - input[j], zero(T))
    end
    step
end

@inline function inner_multiply(::Type{DiagonalJonesMatrix}, X, Y)
    DiagonalJonesMatrix(X.xx'*Y.xx + X.yx'*Y.yx,
                        X.xy'*Y.xy + X.yy'*Y.yy)
end

@inline inner_multiply(::Type{JonesMatrix}, X, Y) = X'*Y

immutable StefcalStep <: StepFunction end
const stefcalstep = StefcalStep()
call(::StefcalStep, input, measured, model) = step = stefcal_step(input, measured, model)
return_type(::StefcalStep, input) = Vector{eltype(input)}

