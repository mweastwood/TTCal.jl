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

@doc """
Calibrate the given measurement set!
""" ->
function bandpass(ms::Table,
                  sources::Vector{PointSource},
                  criteria::StoppingCriteria;
                  force_imaging_columns::Bool = false,
                  model_already_present::Bool = false,
                  reference_antenna::Int = 1)
    frame = reference_frame(ms)
    u,v,w = uvw(ms)
    ν = freq(ms)
    ant1,ant2 = ants(ms)

    sources = filter(source -> isabovehorizon(frame,source),sources)

    Nfreq = length(ν)
    Nant  = numrows(Table(ms[kw"ANTENNA"]))

    gains = ones(Complex64,Nant,2,Nfreq)
    gain_flags = zeros(Bool,Nant,2,Nfreq)

    data  = ms["DATA"]
    model = model_already_present? ms["MODEL_DATA"] : genvis(frame,sources,u,v,w,ν)
    data_flags = ms["FLAG"]
    row_flags  = ms["FLAG_ROW"]
    for α = 1:length(row_flags)
        if row_flags[α]
            data_flags[:,:,α] = true
        end
    end

    if force_imaging_columns || Tables.checkColumnExists(ms,"MODEL_DATA")
        ms["MODEL_DATA"] = model
    end

    bandpass!(gains,gain_flags,data,model,data_flags,
              ant1,ant2,criteria,reference_antenna)
    gains,gain_flags
end

################################################################################
# Internal Interface

function bandpass!(gains::Array{Complex64,3},
                   gain_flags::Array{Bool,3},
                   data::Array{Complex64,3},
                   model::Array{Complex64,3},
                   data_flags::Array{Bool,3},
                   ant1::Vector{Int32},
                   ant2::Vector{Int32},
                   criteria::StoppingCriteria,
                   reference_antenna::Int)
    # Transpose the data and model arrays to create a better memory access pattern
    data  = permutedims(data, (3,1,2))
    model = permutedims(model,(3,1,2))
    data_flags = permutedims(data_flags,(3,1,2))

    Nfreq = size(gains,3)
    for β = 1:Nfreq
        # Calibrate the X polarization
        bandpass_onechannel!(slice(gains,:,1,β),
                             slice(gain_flags,:,1,β),
                             slice( data,:,1,β),
                             slice(model,:,1,β),
                             slice(data_flags,:,1,β),
                             ant1, ant2,
                             criteria,
                             reference_antenna)
        # Calibrate the Y polarization
        bandpass_onechannel!(slice(gains,:,2,β),
                             slice(gain_flags,:,2,β),
                             slice( data,:,4,β),
                             slice(model,:,4,β),
                             slice(data_flags,:,4,β),
                             ant1, ant2,
                             criteria,
                             reference_antenna)
    end
end

@doc """
Calibrate the complex gains from a single frequency channel using a two step process:

1. If the entire channel is flagged, don't bother calibrating.
   Just flag the solutions and move on.
2. Pack the visibilities into square, Hermitian matrices.
3. Flag the auto-correlations.
4. Get an initial estimate of the complex gains.
5. Iteratively improve that estimate.
6. Fix the phase of the reference antenna.
7. Flag the antennas with no unflagged data.
""" ->
function bandpass_onechannel!(gains, gain_flags,
                              data, model, data_flags,
                              ant1::Vector{Int32},
                              ant2::Vector{Int32},
                              criteria::StoppingCriteria,
                              reference_antenna::Int)
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
    square_data  = square_data  .* !square_flags
    gains[:] = firstguess(square_data,square_model)
    square_model = square_model .* !square_flags

    # Elaborating a little bit more on the application of flags:
    #
    # Notice that firstguess(...) makes use of the matrix
    # square_data./square_model. This means that is enough to
    # zero-out the flagged entries of square_data.
    #
    # However, bandpass_inner(...) separately uses the elements
    # of square_data and square_model. Therefore we need
    # both matrices to have their flagged elements zeroed.

    oldgains = similar(gains)
    workspace = RKWorkspace(oldgains,4)

    # 5. Iteratively improve that estimate.
    iter = 0
    converged = false
    while !converged && iter < criteria.maxiter
        oldgains[:] = gains
        rkstep!(gains,workspace,bandpass_step!,square_data,square_model)
        if vecnorm(gains-oldgains)/vecnorm(oldgains) < criteria.tolerance
            converged = true
        end
        iter += 1
    end

    # 6. Fix the phase of the reference antenna.
    factor = conj(gains[reference_antenna])/abs(gains[reference_antenna])
    for ant = 1:length(gains)
        gains[ant] = gains[ant]*factor
    end

    # 7. Flag the antennas with no unflagged data.
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
function makesquare(input,ant1::Vector{Int32},ant2::Vector{Int32})
    N = length(input)
    M = round(Integer,div(sqrt(1+8N)-1,2))
    output = zeros(eltype(input),M,M)
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
function bandpass_inner!(output,input,V,M)
    N = length(input)
    output[:] = complex(0.)
    normalization = zeros(Float32,N)
    @inbounds for j = 1:N, i = 1:N
        z = conj(input[j])*M[i,j]
        output[i] += conj(z)*V[i,j]
        normalization[i] += abs2(z)
    end
    @inbounds for i = 1:N
        normalization[i] == 0. && continue
        output[i] = output[i]/normalization[i] - input[i]
    end
    nothing
end
const bandpass_step! = RKInnerStep{:bandpass_inner!}()

function χ2(data,model,flags,gains)
    output = 0.0
    for i = 1:length(gains), j = 1:i-1
        output += (1-flags[i,j])*abs2(data[i,j] - gains[i]*conj(gains[j])*model[i,j])
    end
    output
end

