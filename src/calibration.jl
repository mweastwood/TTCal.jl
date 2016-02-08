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
            jones::Matrix{$jones}
            flags::Matrix{Bool}
        end
        function $name(Nant, Nfreq)
            $name(ones($jones, Nant, Nfreq), zeros(Bool, Nant, Nfreq))
        end
        Base.similar(cal::$name) = $name(Nant(cal), Nfreq(cal))
    end |> esc
end

@generate_calibration GainCalibration DiagonalJonesMatrix
@generate_calibration PolarizationCalibration JonesMatrix

Nant( cal::Calibration) = size(cal.jones, 1)
Nfreq(cal::Calibration) = size(cal.jones, 2)

doc"""
    invert(calibration::Calibration)

Returns the inverse of the given calibration.
The Jones matrix $J$ of each antenna is set to $J^{-1}$.
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

write(filename,calibration::Calibration) = JLD.save(File(format"JLD",filename),"cal",calibration)
read(filename) = JLD.load(filename,"cal")

function write_for_python(filename, calibration::GainCalibration)
    gains = zeros(Complex128,2, Nant(calibration), Nfreq(calibration))
    flags = zeros(      Bool,   Nant(calibration), Nfreq(calibration))
    for β = 1:Nfreq(calibration), ant = 1:Nant(calibration)
        gains[1,ant,β] = calibration.jones[ant,β].xx
        gains[2,ant,β] = calibration.jones[ant,β].yy
        flags[  ant,β] = calibration.flags[ant,β]
    end
    npzwrite(filename, Dict("gains" => gains, "flags" => flags))
end

# corrupt / applycal

"""
    corrupt!(visibilities::Visibilities, meta::Metadata, calibration::Calibration)

Corrupt the visibilities as if they were observed with the given calibration.
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


"""
    applycal!(visibilities::Visibilities, meta::Metadata, calibration::Calibration)

Apply the calibration to the given visibilities.
"""
function applycal!(visibilities::Visibilities, meta::Metadata, calibration::Calibration)
    inverse_cal = invert(calibration)
    corrupt!(visibilities, meta, inverse_cal)
    visibilities
end

# gaincal / polcal

#=
const argument_docs = """
**Arguments:**

* `ms` - the measurement set from which to derive the calibration
* `sources` - the list of points sources to use as the sky model
* `beam` - the beam model

**Keyword Arguments:**

* `maxiter` - the maximum number of Runge-Kutta steps to take on each
    frequency channel
* `tolerance` - the relative tolerance to use while checking to see if
    more iterations are required
* `flag` - if set to true, attempt to identify and flag slowly converging
    calibration solutions
* `minuvw` - the minimum baseline length (measured in wavelengths) to be
    used during the calibration procedure
* `reference_antenna` - a string containing the antenna number and polarization
    whose phase will be chosen to be zero (eg. "14y" or "62x")
* `force_imaging_columns` - if this is set to true, the MODEL_DATA column
    will be created and populated with model visibilities even if it
    doesn't already exist
"""

"""
    gaincal(ms::MeasurementSet, sources::Vector{Source}, beam::BeamModel;
            maxiter = 20, tolerance = 1e-3, flag = false, minuvw = 0.0,
            reference_antenna = "1x", force_imaging_columns = false)

Solve for the interferometer's electronic gains.

$argument_docs
"""
function gaincal(ms::MeasurementSet, sources::Vector{Source}, beam::BeamModel;
                 maxiter::Int = 20, tolerance::Float64 = 1e-3, flag::Bool = false,
                 minuvw::Float64 = 0.0, reference_antenna::ASCIIString = "1x",
                 force_imaging_columns::Bool = false)
    sources = abovehorizon(ms.frame, sources)
    calibration = GainCalibration(ms.Nant, ms.Nfreq)
    data  = get_data(ms)
    model = genvis(ms, sources, beam)
    flags = get_flags(ms)
    flag_short_baselines!(flags, minuvw, ms.u, ms.v, ms.w, ms.ν)
    set_model_data!(ms, model, force_imaging_columns)
    solve!(calibration, data,model, flags,
           ms.ant1, ms.ant2, maxiter, tolerance, flag)
    fixphase!(calibration, reference_antenna)
    calibration
end

"""
    polcal(ms::MeasurementSet, sources::Vector{Source}, beam::BeamModel;
           maxiter = 20, tolerance = 1e-3, flag = false, minuvw = 0.0,
           reference_antenna = "1x", force_imaging_columns = false)

Solve for the polarization properties of the interferometer.

$argument_docs
"""
function polcal(ms::MeasurementSet, sources::Vector{Source}, beam::BeamModel;
                maxiter::Int = 20, tolerance::Float64 = 1e-3, flag::Bool = false,
                minuvw::Float64 = 0.0, reference_antenna::ASCIIString = "1x",
                force_imaging_columns::Bool = false)
    sources = abovehorizon(ms.frame, sources)
    calibration = PolarizationCalibration(ms.Nant, ms.Nfreq)
    data  = get_corrected_data(ms)
    model = genvis(ms, sources, beam)
    flags = get_flags(ms)
    flag_short_baselines!(flags, minuvw, ms.u, ms.v, ms.w, ms.ν)
    set_model_data!(ms, model)
    solve!(calibration, data, model, flags,
           ms.ant1, ms.ant2, maxiter, tolerance, flag)
    fixphase!(calibration, reference_antenna)
    calibration
end

function solve!(cal::Calibration, data, model, flags, ant1, ant2,
                maxiter, tolerance, switch; quiet::Bool = false)
    quiet || (p = Progress(Nfreq(cal)))
    for β = 1:Nfreq(cal)
        solve_onechannel!(slice(cal.jones,:,β),
                          slice(cal.flags,:,β),
                          slice(data, :,β,:),
                          slice(model,:,β,:),
                          slice(flags,:,β,:),
                          ant1, ant2, maxiter, tolerance, switch)
        quiet || next!(p)
    end
    if !quiet
        # Print a summary of the flags
        print("Flagging summary: ")
        percentage = [sum(cal.flags[ant,:]) / Nfreq(cal) for ant = 1:Nant(cal)]
        idx = percentage .> 0.25
        antennas   = (1:Nant(cal))[idx]
        percentage = percentage[idx]
        if length(antennas) > 0
            for i = 1:length(antennas)
                color = :white
                percentage[i]  ≥ 0.5 && (color = :red)
                percentage[i] == 1.0 && (color = :blue)
                print_with_color(color, string(antennas[i]))
                i == length(antennas) || print(", ")
            end
            print("\n")
            print("          Legend: ")
            print_with_color(:white, ">25% flagged"); print(", ")
            print_with_color(  :red, ">50% flagged"); print(", ")
            print_with_color( :blue, "100% flagged"); print("\n")
        else
            println("all antennas have <25% flags")
        end
    end
end

function solve_onechannel!(jones, jones_flags, data, model, data_flags,
                           ant1, ant2, maxiter, tolerance, switch)
    # If the entire channel is flagged, don't bother calibrating.
    all(data_flags) && (jones_flags[:] = true; return)

    N = length(jones)
    T = eltype(jones)
    square_data  = makesquare( data, data_flags, ant1, ant2)
    square_model = makesquare(model, data_flags, ant1, ant2)
    best_jones   = copy(jones)

    converged = iterate(stefcal, RK4, maxiter, tolerance, switch,
                        best_jones, square_data, square_model)

    jones[:] = best_jones
    jones_flags[:] = solution_flags(best_jones, square_data, square_model, converged)
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
    output = zeros(JonesMatrix, Nant, Nant)
    for α = 1:Nbase
        (flags[1,α] || flags[2,α] || flags[3,α] || flags[4,α]) && continue
        ant1[α] == ant2[α] && continue

        V = JonesMatrix(data[1,α], data[2,α],
                        data[3,α], data[4,α])
        output[ant1[α],ant2[α]] = V
        output[ant2[α],ant1[α]] = output[ant1[α],ant2[α]]'
    end
    output
end

function solution_flags(jones, data, model, converged)
    Nant  = length(jones)
    flags = fill(false, Nant)

    # Flag everything if the solution did not converge.
    if !converged
        flags[:] = true
        return flags
    end

    # Find flagged antennas by looking for columns of zeros.
    for j = 1:Nant
        isflagged = true
        for i = 1:Nant
            if data[i,j] != zero(JonesMatrix)
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

# step functions

doc"""
    stefcal_step(input, data, model) -> step

Given the `data` and `model` visibilities, and the current
guess for the Jones matrices, solve for `step` such
that the new value of the Jones matrices is `input+step`.

The update step is defined such that the new value of the
Jones matrices minimizes

\\[
    \sum_{i,j}\|V_{i,j} - G_i M_{i,j} G_{j,\rm new}^*\|^2,
\\]

where $i$ and $j$ label the antennas, $V$ labels the measured
visibilities, $M$ labels the model visibilities, and $G$
labels the Jones matrices.

*References:*

* Michell, D. et al. 2008, JSTSP, 2, 5.
* Salvini, S. & Wijnholds, S. 2014, A&A, 571, 97.
"""
function stefcal_step{T}(input::AbstractVector{T}, data, model)
    Nant = length(input)
    step = similar(input)
    @inbounds for j = 1:Nant
        numerator   = zero(T)
        denominator = zero(T)
        for i = 1:Nant
            GM = input[i]*model[i,j]
            V  = data[i,j]
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

immutable StefCalStep <: StepFunction end
const stefcal = StefCalStep()
call(::StefCalStep, input, data, model) = step = stefcal_step(input, data, model)

function check!(::StefCalStep, input, data, model)
    Nant = length(input)

    # There are two types of residuals that we can compute:
    #
    # 1. |Vᵢⱼ - GᵢMᵢⱼGⱼ|, and
    # 2. |Gᵢ⁻¹VᵢⱼGⱼ⁻¹ - Mᵢⱼ|.
    #
    # The calibration algorithm used by TTCal attempts to minimize
    # the first kind of residual. In my experience TTCal is very good
    # at finding a solution where the first residual is low even when
    # there exists a number of antennas that are feeding garbage into
    # the algorithm and should be flagged.
    #
    # For example an antenna that is unplugged and is simply picking
    # up noise can be made to minimize the first residual by dialing
    # down the gain amplitudes. However this will make the second
    # residual large.
    #
    # Therefore, let's look for antennas where the second kind of
    # residual is exploding.

    residual = zeros(Nant, Nant)
    for j = 1:Nant, i = j+1:Nant
        data[i,j] == zero(JonesMatrix) && continue
        residual[i,j] = norm(input[i]\data[i,j]/input[j]' - model[i,j])
        residual[j,i] = residual[i,j]
    end

    δ = zeros(Nant)
    for i = 1:Nant
        δ[i] = median(slice(residual, :, i))
    end

    worst_antenna = indmax(δ)
    if δ[worst_antenna] > 10median(δ)
        data[worst_antenna,:] = zero(JonesMatrix)
        data[:,worst_antenna] = zero(JonesMatrix)
        model[worst_antenna,:] = zero(JonesMatrix)
        model[:,worst_antenna] = zero(JonesMatrix)
    end

    nothing
end
=#

