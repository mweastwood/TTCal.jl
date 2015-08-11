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

immutable GainCalibration <: Calibration
    gains::Array{Complex64,3}
    flags::Array{Bool,3}
end

function GainCalibration(Nant,Nfreq)
    gains = ones(Complex64,Nant,Nfreq,2)
    flags = zeros(Bool,Nant,Nfreq,2)
    GainCalibration(gains,flags)
end

Nant(g::GainCalibration) = size(g.gains,1)
Nfreq(g::GainCalibration) = size(g.gains,2)

function invert!(g::GainCalibration)
    for i in eachindex(g.gains)
        g.gains[i] = inv(g.gains[i])
    end
end

function invert(g::GainCalibration)
    out = deepcopy(g)
    invert!(out)
    out
end

"""
Set the phase of the reference antenna to zero.
"""
function fixphase!(g::GainCalibration,
                   reference_antenna)
    for pol = 1:2, β = 1:Nfreq(g)
        factor = (conj(g.gains[reference_antenna,β,pol])
                    / abs(g.gains[reference_antenna,β,pol]))
        for ant = 1:Nant(g)
            g.gains[ant,β,pol] = g.gains[ant,β,pol]*factor
        end
    end
end

@doc """
Calibrate the given measurement set!
""" ->
function bandpass(ms::Table,
                  sources::Vector{PointSource};
                  maxiter::Int = 20,
                  tolerance::Float64 = 1e-5,
                  force_imaging_columns::Bool = false,
                  reference_antenna::Int = 1)
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
    calibration = GainCalibration(Nant,Nfreq)

    data  = ms["DATA"]
    model = genvis(frame,phase_dir,sources,u,v,w,ν)
    flags = MeasurementSets.flags(ms)

    if force_imaging_columns || Tables.checkColumnExists(ms,"MODEL_DATA")
        ms["MODEL_DATA"] = model
    end

    bandpass_internal!(calibration,data,model,flags,
                       ant1,ant2,maxiter,tolerance,
                       reference_antenna)
    calibration
end

################################################################################
# Internal Interface

function bandpass_internal!(calibration,data,model,flags,
                            ant1,ant2,maxiter,tolerance,
                            reference_antenna)
    # Re-order the data to make iterating over it faster
    data  = permutedims(data, (3,2,1))
    flags = permutedims(flags,(3,2,1))
    model = permutedims(model,(3,2,1))

    for β = 1:Nfreq(calibration)
        # Calibrate the X polarization
        bandpass_onechannel!(slice(calibration.gains,:,β,1),
                             slice(calibration.flags,:,β,1),
                             slice(data, :,β,1),
                             slice(model,:,β,1),
                             slice(flags,:,β,1),
                             ant1, ant2,
                             maxiter, tolerance)
        # Calibrate the Y polarization
        bandpass_onechannel!(slice(calibration.gains,:,β,2),
                             slice(calibration.flags,:,β,2),
                             slice(data, :,β,4),
                             slice(model,:,β,4),
                             slice(flags,:,β,4),
                             ant1, ant2,
                             maxiter, tolerance)
    end

    fixphase!(calibration,reference_antenna)
end

@doc """
Calibrate the complex gains from a single frequency channel using a two step process:

1. If the entire channel is flagged, don't bother calibrating.
   Just flag the solutions and move on.
2. Pack the visibilities into square, Hermitian matrices.
3. Flag the auto-correlations.
4. Get an initial estimate of the complex gains.
5. Iteratively improve that estimate.
6. Flag the entire channel if the iteration didn't converge.
8. Flag the antennas with no unflagged data.
""" ->
function bandpass_onechannel!(gains, gain_flags,
                              data, model, data_flags,
                              ant1, ant2,
                              maxiter, tolerance)
    # 1. If the entire channel is flagged, don't bother calibrating.
    #    Just flag the solutions and move on.
    if all(data_flags)
        gain_flags[:] = true
        return
    end

    # 2. Pack the visibilities into square, Hermitian matrices.
    square_data  = makesquare(data, ant1,ant2)
    square_model = makesquare(model,ant1,ant2)
    square_flags = makesquare(data_flags,ant1,ant2)

    # 3. Flag the auto-correlations.
    for i = 1:size(square_flags,1)
        square_flags[i,i] = true
    end

    # 4. Get an initial estimate of the complex gains.
    # (note that the model data is flagged after taking this
    # initial estimate to prevent a divide-by-zero error)
    square_data[square_flags] = 0
    gains[:] = firstguess(square_data,square_model)
    square_model[square_flags] = 0

    # Elaborating a little bit more on the application of flags:
    #
    # Notice that firstguess(...) makes use of the matrix
    # square_data./square_model. This means that is enough to
    # zero-out the flagged entries of square_data.
    #
    # However, bandpass_inner(...) separately uses the elements
    # of square_data and square_model. Therefore we need
    # both matrices to have their flagged elements zeroed.

    # 5. Iteratively improve that estimate.
    converged = @iterate(BandpassStep(),RK4,maxiter,tolerance,
                         gains,square_data,square_model)

    # 6. Flag the entire channel if the iteration didn't converge.
    if !converged
        gain_flags[:] = true
        return
    end

    # 8. Flag the antennas with no unflagged data.
    bad_ants = squeeze(all(square_flags,1),1)
    gains[bad_ants] = 1.0
    gain_flags[bad_ants] = true

    nothing
end

@doc """
Pack the data into a square Hermitian matrix.
This method assumes that the data set has no baselines
with missing data.
""" ->
function makesquare(input,ant1,ant2)
    N = length(input)
    M = round(Integer,div(sqrt(1+8N)-1,2))
    output = zeros(eltype(input),M,M)
    makesquare!(output,input,ant1,ant2)
    output
end

function makesquare!(output,input,ant1,ant2)
    N = length(input)
    M = round(Integer,div(sqrt(1+8N)-1,2))
    for α = 1:N
        if ant1[α] == ant2[α]
            output[ant1[α],ant1[α]] = input[α]
        else
            output[ant1[α],ant2[α]] = input[α]
            output[ant2[α],ant1[α]] = conj(input[α])
        end
    end
    output
end

@doc """
Take a rough first guess at the visibilities.

This function works by looking for the principle eigenvector
of the matrix Gij = (Vij/Mij) where V is the matrix of measured
visibilities, and M is the matrix of model visibilities. In the
absence of noise, and if the model is complete, then G = gg',
where g is the vector of complex gains.

The main confounding factor for this method is the presence of
a nonzero noise term on the diagonal (from the autocorrelations).
Most calibration routines neglect the autocorrelations for this
reason, but there is no obvious way to exclude them in this case.
Hence we use this as a first guess to an iterative method that
refines the calibration.
""" ->
function firstguess(V,M)
    G = V./M
    λ,v = eigs(G,nev=1,which=:LM)
    w = sqrt(λ[1])
    v.*w
end

@doc """
This function defines the step for an iterative method of
solving for the complex gains. It is based on the work of
Stef Salvini.

This method avoids the need to construct and invert huge
Jacobian matrices by ignoring the derivative of the residuals
altogether. Instead an incredibly simple update step is
defined. Convergence is then accelerated by taking clever
linear combinations of these steps (as is done in numerical
ODE solvers).

The update step defined here is g[i] → linear least squares
solution of V = g[i]*M*g', where V and M are vectors containing
the measured and model visibilities where antenna i is the
first antenna.
""" ->
function bandpass_step(g,V,M)
    N = length(g)
    output = zeros(Complex64,N)
    normalization = zeros(Float32,N)
    @inbounds for j = 1:N, i = 1:N
        z = conj(g[j])*M[i,j]
        output[i] += conj(z)*V[i,j]
        normalization[i] += abs2(z)
    end
    @inbounds for i = 1:N
        normalization[i] == 0. && continue
        output[i] = output[i]/normalization[i] - g[i]
    end
    output
end

immutable BandpassStep <: StepFunction end
call(::BandpassStep,g,V,M) = bandpass_step(g,V,M)

