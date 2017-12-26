# Copyright (c) 2015-2017 Michael Eastwood
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

const CalElementTypes = Union{Complex128, DiagonalJonesMatrix, JonesMatrix}

struct Solution{P <: Polarization, T} <: AbstractVector{eltype(T)}
    data :: T
end

function Solution(pol::Type{<:Polarization}, Nant)
    T = vector_type(pol)
    data = T(Nant)
    for antenna = 1:Nant
        data[antenna] = one(eltype(data))
    end
    Solution{pol, T}(data)
end

Nant(cal::Solution) = length(cal.data)
polarization(::Solution{P, T}) where {P, T} = P
Base.eltype( ::Solution{P, T}) where {P, T} = eltype(T)

Base.getindex( cal::Solution, antenna) = cal.data[antenna]
Base.setindex!(cal::Solution, value, antenna) = cal.data[antenna] = value
Base.similar(  cal::Solution) = Solution(polarization(cal), Nant(cal))
Base.length(   cal::Solution) = length(cal.data)
Base.size(     cal::Solution) = size(cal.data)

struct Calibration{P <: Polarization, T}
    data :: Matrix{Solution{P, T}}
end

Nfreq(cal::Calibration) = size(cal.data, 1)
Ntime(cal::Calibration) = size(cal.data, 2)

Base.getindex(cal::Calibration, frequency, time) = cal.data[frequency, time]

function Calibration(metadata; polarization=Full, collapse_frequency=false, collapse_time=false)
    N = collapse_frequency ? 1 : Nfreq(metadata)
    M = collapse_time      ? 1 : Ntime(metadata)
    data = [Solution(polarization, Nant(metadata)) for freq = 1:N, time = 1:M]
    Calibration(data)
end

function calibrate(data::Dataset, model::Dataset;
                   collapse_frequency=false, collapse_time=false,
                   maxiter=50, tolerance=1e-3, minuvw=15.0, quiet=false)
    flag_short_baselines!(data, minuvw)
    calibration = Calibration(data.metadata, polarization=polarization(data),
                              collapse_frequency=collapse_frequency, collapse_time=collapse_time)
    calibrate!(calibration, data, model, maxiter, tolerance, quiet)
    calibration
end

function calibrate!(calibration::Calibration, data::Dataset, model::Dataset,
                    maxiter, tolerance, quiet)
    match_flags!(model, data)
    if Nfreq(calibration) == Ntime(calibration) == 1
        calibrate_onefrequency_onetime!(calibration, data, model, maxiter, tolerance, quiet)
    elseif Nfreq(calibration) == 1
        calibrate_onefrequency!(calibration, data, model, maxiter, tolerance, quiet)
    elseif Ntime(calibration) == 1
        calibrate_onetime!(calibration, data, model, maxiter, tolerance, quiet)
    else
        quiet || (prg = Progress(Ntime(calibration)*Nfreq(calibration)))
        for time = 1:Ntime(calibration), frequency = 1:Nfreq(calibration)
            solve!(calibration[frequency, time].data,
                   data[frequency, time], model[frequency, time],
                   maxiter, tolerance, quiet)
            quiet || next!(prg)
        end
    end
end

function calibrate_onefrequency!(calibration::Calibration, data::Dataset, model::Dataset,
                                 maxiter, tolerance, quiet)
    quiet || (prg = Progress(Ntime(calibration)))
    for time = 1:Ntime(calibration)
        data_slice  = [getindex.( data[:, time], antenna) for antenna = 1:Nant(data)]
        model_slice = [getindex.(model[:, time], antenna) for antenna = 1:Nant(data)]
        solve!(calibration[1, time].data, data_slice, model_slice, maxiter, tolerance, quiet)
        quiet || next!(prg)
    end
end

function calibrate_onetime!(calibration::Calibration, data::Dataset, model::Dataset,
                            maxiter, tolerance, quiet)
    quiet || (prg = Progress(Nfreq(calibration)))
    for frequency = 1:Nfreq(calibration)
        data_slice  = [getindex.( data[frequency, :], antenna) for antenna = 1:Nant(data)]
        model_slice = [getindex.(model[frequency, :], antenna) for antenna = 1:Nant(data)]
        solve!(calibration[frequency, 1].data, data_slice, model_slice, maxiter, tolerance, quiet)
        quiet || next!(prg)
    end
end

function calibrate_onefrequency_onetime!(calibration::Calibration, data::Dataset, model::Dataset,
                                         maxiter, tolerance, quiet)
    data_slice  = [getindex.( data[:], antenna) for antenna = 1:Nant(data)]
    model_slice = [getindex.(model[:], antenna) for antenna = 1:Nant(data)]
    solve!(calibration[1, 1].data, data_slice, model_slice, maxiter, tolerance, quiet)
end

function solve!(gains, measured_visibilities, model_visibilities, maxiter, tolerance, quiet)
    workspace = RKWorkspace(gains, 4)
    function step!(output, input)
        stefcal_step!(output, input, measured_visibilities, model_visibilities)
    end
    iter = 0
    converged = false
    while !converged && iter < maxiter
        newgains = rk4!(workspace, step!, gains)
        relative_change = norm(gains - newgains)/norm(gains)
        if relative_change < tolerance
            converged = true
        end
        gains[:] = newgains
        iter += 1
    end
    if !quiet && !converged
        warn("calibration did not converge")
    end
    gains
end







#"""
#    GainCalibration
#
#*Description*
#
#This type represents the gain calibration of an interferometer.
#
#Each antenna and frequency channel receives a diagonal Jones matrix.
#The diagonal terms of the Jones matrix are the complex gains of the `x`
#and `y` polarization respectively. The off-diagonal terms represent the
#polarization leakage from `x` to `y` and `y` to `x`. These off-diagonal
#terms are assumed to be zero.
#
#The `gaincal` routine is used to solve for the interferometer's
#gain calibration.
#
#*Fields*
#
#* `jones` - an array of diagonal Jones matrices (one per antenna and frequency channel)
#* `flags` - a corresponding list of flags
#"""
#@generate_calibration GainCalibration DiagonalJonesMatrix
#
#"""
#    PolarizationCalibration
#
#*Description*
#
#This type represents the polarization calibration of an interferometer.
#
#Each antenna and frequency channel receives a full Jones matrix.
#The diagonal terms of the Jones matrix are the complex gains of the `x`
#and `y` polarization respectively. The off-diagonal terms represent the
#polarization leakage from `x` to `y` and `y` to `x`. All of these
#terms are included in this calibration.
#
#The `polcal` routine is used to solve for the interferometer's
#polarization calibration.
#
#*Fields*
#
#* `jones` - an array of Jones matrices (one per antenna and frequency channel)
#* `flags` - a corresponding list of flags
#"""
#@generate_calibration PolarizationCalibration JonesMatrix
#
#Nant( cal::Calibration) = size(cal.jones, 1)
#Nfreq(cal::Calibration) = size(cal.jones, 2)
#
#write(filename, calibration::Calibration) = JLD.save(File(format"JLD", filename), "cal", calibration)
#read(filename) = JLD.load(filename, "cal")
#
## The following functions are disabled until NPZ.jl is fixed on Julia v0.5. To re-enable:
## * uncomment these lines
## * add `using NPZ` to `src/TTCal.jl`
## * add `NPZ` to `REQUIRE`
## * uncomment `using NPZ` in `test/runtests.jl`
## * uncomment the corresponding tests in `test/calibration.jl`
##
##function write_for_python(filename, calibration::GainCalibration)
##    gains = zeros(Complex128, 2, Nant(calibration), Nfreq(calibration))
##    flags = zeros(      Bool,    Nant(calibration), Nfreq(calibration))
##    for β = 1:Nfreq(calibration), ant = 1:Nant(calibration)
##        gains[1,ant,β] = calibration.jones[ant,β].xx
##        gains[2,ant,β] = calibration.jones[ant,β].yy
##        flags[  ant,β] = calibration.flags[ant,β]
##    end
##    npzwrite(filename, Dict("gains" => gains, "flags" => flags))
##end
##
##function write_for_python(filename, calibration::PolarizationCalibration)
##    gains = zeros(Complex128, 4, Nant(calibration), Nfreq(calibration))
##    flags = zeros(      Bool,    Nant(calibration), Nfreq(calibration))
##    for β = 1:Nfreq(calibration), ant = 1:Nant(calibration)
##        gains[1,ant,β] = calibration.jones[ant,β].xx
##        gains[2,ant,β] = calibration.jones[ant,β].xy
##        gains[3,ant,β] = calibration.jones[ant,β].yx
##        gains[4,ant,β] = calibration.jones[ant,β].yy
##        flags[  ant,β] = calibration.flags[ant,β]
##    end
##    npzwrite(filename, Dict("gains" => gains, "flags" => flags))
##end
#
## gaincal / polcal
#
#"""
#    gaincal(visibilities, metadata, beam, sources)
#    gaincal(visibilities, metadata, model_visibilities)
#
#*Description*
#
#Solve for the interferometer's electronic gains.
#
#*Arguments*
#
#* `visibilities` - the visibilities measured by the interferometer
#* `metadata` - the metadata describing the interferometer
#* `beam` - the primary beam model
#* `sources` - a list of sources comprising the sky model
#* `model_visibilities` - alternatively the sky model visibilities can be provided
#
#*Keyword Arguments*
#
#* `maxiter` - the maximum number of iterations to take on each frequency channel (defaults to `20`)
#* `tolerance` - the relative tolerance used to test for convergence (defaults to `1e-3`)
#* `quiet` - suppresses printing if set to `true` (defaults to `false`)
#"""
#function gaincal{S<:Source}(visibilities::Visibilities, meta::Metadata, beam::BeamModel, sources::Vector{S};
#                            maxiter = 20, tolerance = 1e-3, quiet = false)
#    frame = reference_frame(meta)
#    sources = abovehorizon(frame, sources)
#    model = genvis(meta, beam, sources)
#    gaincal(visibilities, meta, model, maxiter=maxiter, tolerance=tolerance, quiet=quiet)
#end
#
#function gaincal(visibilities::Visibilities, meta::Metadata, model::Visibilities;
#                 maxiter = 20, tolerance = 1e-3, quiet = false)
#    calibration = GainCalibration(Nant(meta), Nfreq(meta))
#    solve!(calibration, visibilities, model, meta, maxiter, tolerance, quiet)
#    calibration
#end
#
#"""
#    polcal(visibilities, metadata, beam, sources)
#    polcal(visibilities, metadata, model_visibilities)
#
#*Description*
#
#Solve for the polarization properties of the interferometer.
#
#*Arguments*
#
#* `visibilities` - the visibilities measured by the interferometer
#* `metadata` - the metadata describing the interferometer
#* `sources` - a list of sources comprising the sky model
#* `beam` - the primary beam model
#* `model_visibilities` - alternatively the sky model visibilities can be provided
#
#*Keyword Arguments*
#
#* `maxiter` - the maximum number of iterations to take on each frequency channel (defaults to `20`)
#* `tolerance` - the relative tolerance used to test for convergence (defaults to `1e-3`)
#* `quiet` - suppresses printing if set to `true` (defaults to `false`)
#"""
#function polcal{S<:Source}(visibilities::Visibilities, meta::Metadata, beam::BeamModel, sources::Vector{S};
#                           maxiter::Int = 20, tolerance::Float64 = 1e-3, quiet = false)
#    frame = reference_frame(meta)
#    sources = abovehorizon(frame, sources)
#    model = genvis(meta, beam, sources)
#    polcal(visibilities, meta, model, maxiter=maxiter, tolerance=tolerance, quiet=quiet)
#end
#
#function polcal(visibilities::Visibilities, meta::Metadata, model::Visibilities;
#                maxiter = 20, tolerance = 1e-3, quiet = false)
#    calibration = PolarizationCalibration(Nant(meta), Nfreq(meta))
#    solve!(calibration, visibilities, model, meta, maxiter, tolerance, quiet)
#    calibration
#end
#
#function solve!(calibration, measured_visibilities, model_visibilities, metadata, maxiter, tolerance, quiet)
#    println(2)
#    square_measured, square_model = makesquare(measured_visibilities, model_visibilities, metadata)
#    println(2)
#    quiet || (p = Progress(Nfreq(metadata), "Calibrating: "))
#    for β = 1:Nfreq(metadata)
#        @show β
#        solve_onechannel!(view(calibration.jones, :, β),
#                          view(calibration.flags, :, β),
#                          view(square_measured, :, :, β),
#                          view(square_model,    :, :, β),
#                          maxiter, tolerance)
#        quiet || next!(p)
#    end
#    if !quiet
#        # Print a summary of the flags
#        flagged_channels = sum(all(calibration.flags, 1))
#        flagged_antennas = sum(all(calibration.flags, 2))
#        percentage = 100 * sum(calibration.flags) / length(calibration.flags)
#        @printf("(%d antennas flagged, %d channels flagged, %0.2f percent total)\n",
#                flagged_antennas, flagged_channels, percentage)
#    end
#end
#
#"Solve one channel at a time."
#function solve_onechannel!(jones, flags, measured, model, maxiter, tolerance)
#    converged = iterate(RK4(stefcal_step), maxiter, tolerance, jones, measured, model)
#    flag_solution!(jones, flags, measured, model, converged)
#    jones
#end
#
#"Solve with one solution for all channels."
#function solve_allchannels!(calibration, measured, model, metadata, maxiter, tolerance)
#    G = view(calibration.jones, :, 1)
#    F = view(calibration.flags, :, 1)
#    V, M = makesquare(measured, model, metadata)
#    converged = iterate(RK4(stefcal_step), maxiter, tolerance, G, V, M)
#    flag_solution!(G, F, V, M, converged)
#    calibration
#end
#
#doc"""
#Pack the visibilities into a square Hermitian matrices such that
#the visibilities are ordered as follows:
#
#$$\begin{pmatrix}
#V_{11} & V_{12} & V_{13} &        & \\\\
#V_{21} & V_{22} & V_{23} &        & \\\\
#V_{31} & V_{32} & V_{33} &        & \\\\
#       &        &        & \ddots & \\\\
#\end{pmatrix}$$
#
#Auto-correlations and flagged cross-correlations are set to zero.
#"""
#function makesquare(measured, model, meta)
#    output_measured = zeros(JonesMatrix, Nant(meta), Nant(meta), Nfreq(meta))
#    output_model    = zeros(JonesMatrix, Nant(meta), Nant(meta), Nfreq(meta))
#    for β = 1:Nfreq(meta), α = 1:Nbase(meta)
#        measured.flags[α,β] && continue
#        antenna1 = meta.baselines[α].antenna1
#        antenna2 = meta.baselines[α].antenna2
#        antenna1 == antenna2 && continue
#        output_measured[antenna1,antenna2,β] = measured.data[α,β]
#        output_measured[antenna2,antenna1,β] = measured.data[α,β]'
#        output_model[antenna1,antenna2,β] = model.data[α,β]
#        output_model[antenna2,antenna1,β] = model.data[α,β]'
#    end
#    output_measured, output_model
#end
#
#"""
#Inspect the calibration solution and apply flags as necessary.
#"""
#function flag_solution!(jones, flags, measured, model, converged)
#    Nant  = length(jones)
#
#    # Flag everything if the solution did not converge.
#    if !converged
#        flags[:] = true
#        return flags
#    end
#
#    # Find flagged antennas by looking for columns of zeros.
#    for j = 1:Nant
#        isflagged = true
#        for i = 1:Nant
#            if measured[i,j] != zero(JonesMatrix)
#                isflagged = false
#                break
#            end
#        end
#        if isflagged
#            flags[j] = true
#        end
#    end
#
#    flags
#end
#
#"""
#    fixphase!(calibration, reference_antenna)
#
#**Description**
#
#Set the phase of the reference antenna and polarization to zero.
#
#*Arguments*
#
#* `calibration` - the calibration that will have its phase adjusted
#* `reference_antenna` - a string containing the antenna number and polarization
#    whose phase will be chosen to be zero (eg. "14y" or "62x")
#"""
#function fixphase!(cal::Calibration, reference_antenna)
#    regex = r"(\d+)(x|y)"
#    m = match(regex,reference_antenna)
#    refant = parse(Int,m.captures[1])
#    refpol = m.captures[2] == "x"? 1 : 2
#    for β = 1:Nfreq(cal)
#        ref = refpol == 1? cal.jones[refant,β].xx : cal.jones[refant,β].yy
#        factor = conj(ref) / abs(ref)
#        for ant = 1:Nant(cal)
#            cal.jones[ant,β] = cal.jones[ant,β]*factor
#        end
#    end
#    cal
#end

# step functions

function stefcal_step(gains, measured_visibilities, model_visibilities)
    step = similar(gains)
    stefcal_step!(step, gains, measured_visibilities, model_visibilities)
    step
end

function stefcal_step!(steps, gains, measured_visibilities, model_visibilities)
    N = length(gains)
    for antenna = 1:N
        steps[antenna] = _stefcal_step(gains,
                                       measured_visibilities[antenna],
                                       model_visibilities[antenna],
                                       antenna)
    end
    steps
end

function _stefcal_step(gains, measured_visibilities, model_visibilities, antenna)
    numerator, denominator = numerator_denominator(gains,
                                                   measured_visibilities,
                                                   model_visibilities,
                                                   antenna)
    ok = abs(det(denominator)) > eps(Float64)
    if !ok
        return zero(typeof(numerator))
    else
        newgain = (denominator\numerator)'
        oldgain = gains[antenna]
        return newgain - oldgain
    end
end

"handle the case where we are given multiple frequencies or times to solve jointly"
function numerator_denominator(gains                 :: AbstractVector{T},
                               measured_visibilities :: AbstractArray{<:AbstractVector},
                               model_visibilities    :: AbstractArray{<:AbstractVector},
                               antenna               :: Int) where {T<:CalElementTypes}
    numerator   = zero(T)
    denominator = zero(T)
    for idx in eachindex(measured_visibilities, model_visibilities)
        V = measured_visibilities[idx]
        M =    model_visibilities[idx]
        num, den = numerator_denominator(gains, V, M, antenna)
        numerator   += num
        denominator += den
    end
    numerator, denominator
end

"handle the case where we are given a single frequency and time to solve"
function numerator_denominator(gains                 :: AbstractVector{T},
                               measured_visibilities :: AbstractVector{T},
                               model_visibilities    :: AbstractVector{T},
                               antenna               :: Int) where {T<:CalElementTypes}
    N = length(gains)
    numerator   = zero(T)
    denominator = zero(T)
    @simd for idx = 1:N
        @inbounds G = gains[idx]
        @inbounds V = measured_visibilities[idx]
        @inbounds M = model_visibilities[idx]
        GM  = G*M
        GM′ = GM'
        numerator   += GM′* V
        denominator += GM′*GM
    end
    numerator, denominator
end

#doc"""
#    stefcal_step(input, measured, model)
#
#**Description**
#
#Given the `measured` and `model` visibilities, and the current
#guess for the Jones matrices, solve for `step` such
#that the new value of the Jones matrices is `input+step`.
#
#The update step is defined such that the new value of the
#Jones matrices minimizes
#
#$$\sum_{i,j}\|V_{i,j} - J_i M_{i,j} J_{j,\rm new}^*\|^2,$$
#
#where $i$ and $j$ label the antennas, $V$ labels the measured
#visibilities, $M$ labels the model visibilities, and $J$
#labels the Jones matrices.
#
#*References*
#
#* Michell, D. et al. 2008, JSTSP, 2, 5.
#* Salvini, S. & Wijnholds, S. 2014, A&A, 571, 97.
#"""
#function stefcal_step{T}(input::AbstractVector{T}, measured::Matrix, model::Matrix)
#    Nant = length(input) # number of antennas
#    step = similar(input)
#    @inbounds for j = 1:Nant
#        numerator   = zero(T)
#        denominator = zero(T)
#        for i = 1:Nant
#            GM = input[i]*model[i,j]
#            V  = measured[i,j]
#            numerator   += inner_multiply(T, GM, V)
#            denominator += inner_multiply(T, GM, GM)
#        end
#        ok = abs(det(denominator)) > eps(Float64)
#        step[j] = ifelse(ok, (denominator\numerator)' - input[j], zero(T))
#    end
#    step
#end
#
#function stefcal_step{T}(input::AbstractVector{T}, measured, model)
#    # This version is called if `measured` and `model` have 3 dimensions
#    # indicating that we want to use multiple integrations or multiple
#    # frequency channels to solve for the calibration.
#    Nant = length(input)     # number of antennas
#    Nint = size(measured, 3) # number of integrations
#    @show "stefcal", Nant, Nint
#    step = similar(input)
#    @inbounds for j = 1:Nant
#        numerator   = zero(T)
#        denominator = zero(T)
#        for t = 1:Nint, i = 1:Nant
#            GM = input[i]*model[i,j,t]
#            V  = measured[i,j,t]
#            numerator   += inner_multiply(T, GM, V)
#            denominator += inner_multiply(T, GM, GM)
#        end
#        ok = abs(det(denominator)) > eps(Float64)
#        step[j] = ifelse(ok, (denominator\numerator)' - input[j], zero(T))
#    end
#    step
#end
#
#@inline function inner_multiply(::Type{DiagonalJonesMatrix}, X, Y)
#    DiagonalJonesMatrix(X.xx'*Y.xx + X.yx'*Y.yx,
#                        X.xy'*Y.xy + X.yy'*Y.yy)
#end
#
#@inline inner_multiply(::Type{JonesMatrix}, X, Y) = X'*Y

# corrupt / applycal

#doc"""
#    corrupt!(visibilities, metadata, calibration)
#
#*Description*
#
#Corrupt the visibilities as if they were observed with the given calibration.
#That is if we have the model visibility $V_{i,j}$ on baseline $i,j$, and the
#corresponding Jones matrices $J_i$ and $J_j$ for antennas $i$ and $j$ compute
#
#$$V_{i,j} \rightarrow J_i V_{i,j} J_j^*$$
#
#*Arguments*
#
#* `visibilities` - the list of visibilities to corrupt
#* `metadata` - the metadata describing the interferometer
#* `calibration` - the calibration in question
#"""

function corrupt(dataset::Dataset, calibration::Calibration)
    corrupted = deepcopy(dataset)
    corrupt!(corrupted, calibration)
    corrupted
end

function corrupt!(dataset::Dataset, calibration::Calibration)
    if Nfreq(calibration) == Ntime(calibration) == 1
        corrupt_onefrequency_onetime!(dataset, calibration)
    elseif Nfreq(calibration) == 1
        corrupt_onefrequency!(dataset, calibration)
    elseif Ntime(calibration) == 1
        corrupt_onetime!(dataset, calibration)
    else
        for time = 1:Ntime(dataset), frequency = 1:Nfreq(dataset)
            corrupt!(dataset[frequency, time], calibration[frequency, time])
        end
    end
    dataset
end

function corrupt_onefrequency!(dataset::Dataset, calibration::Calibration)
    for time = 1:Ntime(dataset), frequency = 1:Nfreq(dataset)
        corrupt!(dataset[frequency, time], calibration[1, time])
    end
end

function corrupt_onetime!(dataset::Dataset, calibration::Calibration)
    for time = 1:Ntime(dataset), frequency = 1:Nfreq(dataset)
        corrupt!(dataset[frequency, time], calibration[frequency, 1])
    end
end

function corrupt_onefrequency_onetime!(dataset::Dataset, calibration::Calibration)
    for time = 1:Ntime(dataset), frequency = 1:Nfreq(dataset)
        corrupt!(dataset[frequency, time], calibration[1, 1])
    end
end

function corrupt!(visibilities::Visibilities, solution::Solution)
    for antenna1 = 1:Nant(visibilities), antenna2 = antenna1:Nant(visibilities)
        isflagged(visibilities, antenna1, antenna2) && continue
        V  = visibilities[antenna1, antenna2]
        J₁ = solution[antenna1]
        J₂ = solution[antenna2]
        visibilities[antenna1, antenna2] = J₁*V*J₂'
    end
    visibilities
end


#doc"""
#    applycal!(visibilities, metadata, calibration)
#
#*Description*
#
#Apply the calibration to the given visibilities.
#That is if we have the model visibility $V_{i,j}$ on baseline $i,j$, and the
#corresponding Jones matrices $J_i$ and $J_j$ for antennas $i$ and $j$ compute
#
#$$V_{i,j} \rightarrow J_i^{-1} V_{i,j} (J_j^{-1})^*$$
#
#*Arguments*
#
#* `visibilities` - the list of visibilities to corrupt
#* `metadata` - the metadata describing the interferometer
#* `calibration` - the calibration in question
#"""
function applycal(dataset::Dataset, calibration::Calibration)
    calibrated = deepcopy(dataset)
    applycal!(calibrated, calibration)
    calibrated
end

function applycal!(dataset::Dataset, calibration::Calibration)
    corrupt!(dataset, invert(calibration))
end

#doc"""
#    invert(calibration)
#
#Returns the inverse of the given calibration.
#The Jones matrices $J$ of each antenna is set to $J^{-1}$.
#"""
function invert(calibration::Calibration)
    calibration = deepcopy(calibration)
    for time = 1:Ntime(calibration), frequency = 1:Nfreq(calibration)
        solution = calibration[frequency, time]
        for antenna = 1:Nant(solution)
            solution[antenna] = inv(solution[antenna])
        end
    end
    calibration
end

