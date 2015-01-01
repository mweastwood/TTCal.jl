@doc """
This type is a container for the configuration options used in the bandpass
routine.

# Fields
* `maxiter`:
The maximum number of iterations to use in the calibration loop.
* `tol`:
The relative tolerance to use as a convergence criterion. The calibration is
deemed to have converged when the gains change by less than `tol` from the
previous iteration.

# Type Parameters
* `RKn`:
The Runge-Kutta method to use. For example, 4 denotes the RK4 method.
* `doubleprecision`:
A true/false flag deciding whether the complex gains should be 32+32
or 64+64 bit complex numbers.
""" ->
type BandpassOptions{RKn,doubleprecision}
    maxiter::Int
    tol::Float64
end

function BandpassOptions(maxiter::Int,
                         tol::Float64,
                         RKn::Int,
                         doubleprecision::Bool)
    # This is not type stable! Fix!
    BandpassOptions{RKn,doubleprecision}(maxiter,tol)
end

@doc """
This type contains fields for use (and re-use) in the inner calibration loop.
The purpose of this type is to reduce the number of allocations made within
the inner calibration loop by facilitating the re-use large arrays.
""" ->
type BandpassWorkspace{T}
    # Visibilities
    V::Matrix{Complex{T}} # Data
    M::Matrix{Complex{T}} # Model

    # Gains
    g::Vector{Complex{T}}
    newg::Vector{Complex{T}}   # (updated by outer steps)
    trialg::Vector{Complex{T}} # (updated by inner steps)

    # Matrices used in the calibration routine
    normalization::Vector{T}
    k::Matrix{Complex{T}}
end

function getworkspace{RKn,doubleprecision}(interferometer::Interferometer,
                                           bandpass::BandpassOptions{RKn,doubleprecision})
    N = interferometer.Nant - length(interferometer.flaggedantennas) # number of unflagged antennas
    T = doubleprecision? Float64 : Float32

    V  = Array(Complex{T},N,N)
    M  = Array(Complex{T},N,N)
    g  = Array(Complex{T},N)
    newg = Array(Complex{T},N)
    trialg = Array(Complex{T},N)
    normalization = Array(T,N)
    k = Array(Complex{T},N,RKn)

    BandpassWorkspace(V,M,g,newg,trialg,normalization,k)
end

@doc """
Calibrate the given measurement set!
""" ->
function bandpass{_,doubleprecision}(interferometer::Interferometer,
                                     ms::Vector{Table},
                                     sources::Vector{Source},
                                     options::BandpassOptions{_,doubleprecision})
    Nms   = length(ms)
    Nfreq = interferometer.Nfreq
    Nant  = interferometer.Nant
    gains = ones(doubleprecision? Complex128 : Complex64,Nant,2,Nms*Nfreq)
    for i = 1:length(ms)
        bandpass!(sub(gains,:,:,(i-1)*Nfreq+1:i*Nfreq),
                  interferometer,ms[i],sources,options)
    end
    gains
end

function bandpass!(gains,
                   interferometer::Interferometer,
                   ms::Table,
                   sources::Vector{Source},
                   options::BandpassOptions)
    data  = ms["DATA"]
    model = visibilities(ms,sources)
    if Tables.checkColumnExists(ms,"MODEL_DATA")
        ms["MODEL_DATA"] = model
    end
    bandpass!(gains,interferometer,data,model,options)
end

function bandpass!(gains,
                   interferometer::Interferometer,
                   data::Array{Complex64,3},
                   model::Array{Complex64,3},
                   options::BandpassOptions)
    # Transpose the data and model arrays to create a better
    # memory access pattern
    data  = permutedims(data, (3,2,1))
    model = permutedims(model,(3,2,1))

    workspace = getworkspace(interferometer,options)
    # β = frequency channel index
    for β = 1:interferometer.Nfreq
        # TODO: A fully polarized calibration does not deal with the polarizations
        # separately. Fix this!

        # Calibrate the X polarization
        bandpass_onechannel!(slice(gains,:,1,β),
                             slice( data,:,β,1),
                             slice(model,:,β,1),
                             interferometer,
                             options,
                             workspace)
        # Calibrate the Y polarization
        bandpass_onechannel!(slice(gains,:,2,β),
                             slice( data,:,β,4),
                             slice(model,:,β,4),
                             interferometer,
                             options,
                             workspace)
    end
    nothing
end

@doc """
Calibrate the complex gains from a single frequency channel using a two step process:

1. Get an initial estimate of the gains.
2. Iteratively improve that estimate.
""" ->
function bandpass_onechannel!{T}(gains::SubArray{T},
                                 data::SubArray,
                                 model::SubArray,
                                 interferometer::Interferometer,
                                 options::BandpassOptions,
                                 workspace::BandpassWorkspace)
    # Pack the visibilities into square, Hermitian matrices
    makesquare!(workspace.V, data,interferometer.flaggedantennas)
    makesquare!(workspace.M,model,interferometer.flaggedantennas)

    # Take a first guess at the complex gains
    firstguess!(workspace)

    # Refine the estimate of the gains
    iter = 0
    converged = false
    #@show χ2(workspace)
    while !converged && iter < options.maxiter
        outerstep!(workspace,options)
        if vecnorm(workspace.newg-workspace.g)/vecnorm(workspace.g,Inf) < options.tol
            converged = true
        end
        workspace.g = workspace.newg
        iter += 1
        #@show χ2(workspace)
    end
    #@show χ2(workspace)
    #@show iter

    # Output
    idx = 1
    for ant = 1:interferometer.Nant
        ant in interferometer.flaggedantennas && continue
        gains[ant] = workspace.g[idx]
        idx += 1
    end

    # Fix the phase of the reference antenna
    refant = interferometer.refant
    factor = conj(gains[refant])/abs(gains[refant])
    for ant = 1:interferometer.Nant
        ant in interferometer.flaggedantennas && continue
        gains[ant] = gains[ant]*factor
    end
    nothing
end

@doc """
Pack the visibilities into a square matrix under the
assumption that the visibilities are ordered as follows

    11, 12, 13, ..., 22, 23, ..., 33, ...

An optional list of rows/columns to flag can be provided.
""" ->
function makesquare!(output::Matrix,input,flags = Int[])
    N = length(input)
    # Number of rows/columns
    M = round(Integer,div(sqrt(1+8N)-1,2))
    # Number of rows/columns after flagging
    M_ = M - length(flags)
    count   = 1 # current index for input
    i_ = j_ = 1 # current indices for output
    for i = 1:M
        if i in flags
            # Skip a whole row
            count += M+1-i
            continue
        end
        # Do the diagonal
        output[i_,i_] = input[count]
        count += 1
        # Fill out the remainder of the row
        j_ = i_+1
        for j = i+1:M
            if j in flags
                # Skip a column
                count += 1
                continue
            end
            output[i_,j_] =      input[count]
            output[j_,i_] = conj(input[count])
            count += 1
            j_ += 1
        end
        i_ += 1
    end
    nothing
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
function firstguess!(workspace::BandpassWorkspace)
    G = workspace.V./workspace.M
    λ,v  = eigs(G,nev=1,which=:LM)
    workspace.g = v*sqrt(λ)
    nothing
end

@doc """
The Butcher tableau for the 2nd-order Runge-Kutta method.
""" ->
const RK2 = [1/2 0.0
             0.0 1.0]

@doc """
The Butcher tableau for the 3rd-order Runge-Kutta method.
""" ->
const RK3 = [1/2 0.0 0.0;
            -1.0 2.0 0.0;
             1/6 2/3 1/6]

@doc """
The Butcher tableau for the 4th-order Runge-Kutta method.
""" ->
const RK4 = [1/2 0.0 0.0 0.0;
             0.0 1/2 0.0 0.0;
             0.0 0.0 1.0 0.0;
             1/6 1/3 1/3 1/6]

function outerstep!{RKn,_}(w::BandpassWorkspace,
                           options::BandpassOptions{RKn,_})
    # Select the appropriate tableau
    # (the compiler should eliminate the dead branches)
    if     RKn == 2
        tableau = RK2
    elseif RKn == 3
        tableau = RK3
    elseif RKn == 4
        tableau = RK4
    else
        error("Unimplemented Runge-Kutta method.")
    end

    # Take the Runge-Kutta step
    innerstep!(sub(w.k,:,1), w.normalization, w.V, w.M, w.g)
    for iter = 2:RKn
        w.trialg = w.g
        for j = 1:iter-1
            tableau[iter-1,j] == 0.0 && continue
            w.trialg += w.k[:,j] * tableau[iter-1,j]
        end
        innerstep!(sub(w.k,:,iter), w.normalization, w.V, w.M, w.trialg)
    end

    # Output
    w.newg = w.g
    for j = 1:RKn
        w.newg += w.k[:,j] * tableau[end,j]
    end
    nothing
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
function innerstep!(output,normalization,V,M,g)
    N = length(g)
    for i = 1:N
        output[i] = complex(0.)
        normalization[i] = 0.
    end
    for j = 1:N, i = 1:N
        i == j && continue
        z = conj(g[j])*M[i,j]
        output[i] += conj(z)*V[i,j]
        normalization[i] += abs2(z)
    end
    for i = 1:N
        output[i] = output[i]/normalization[i] - g[i]
    end
    nothing
end

@doc """
Calculate the χ² value of the gain calibration.
""" ->
function χ2(data::Matrix,model::Matrix,gains::Vector)
    out = 0.
    for i = 1:length(gains), j = i+1:length(gains)
        out += abs2(data[j,i] - gains[j]*conj(gains[i])*model[j,i])
    end
    out
end

function χ2(workspace::BandpassWorkspace)
    χ2(workspace.V,workspace.M,workspace.g)
end

