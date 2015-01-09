@doc """
Calibrate the given measurement set!
""" ->
function bandpass(interferometer::Interferometer,
                  ms::Vector{Table},
                  sources::Vector{Source},
                  criteria::StoppingCriteria)
    Nms   = length(ms)
    Nfreq = interferometer.Nfreq
    Nant  = interferometer.Nant
    gains = ones(Complex64,Nant,2,Nms*Nfreq)
    for i = 1:length(ms)
        bandpass!(sub(gains,:,:,(i-1)*Nfreq+1:i*Nfreq),
                  interferometer,ms[i],sources,criteria)
    end
    gains
end

function bandpass!(gains,
                   interferometer::Interferometer,
                   ms::Table,
                   sources::Vector{Source},
                   criteria::StoppingCriteria)
    data  = ms["DATA"]
    model = visibilities(ms,sources)
    if Tables.checkColumnExists(ms,"MODEL_DATA")
        ms["MODEL_DATA"] = model
    end
    bandpass!(gains,data,model,interferometer,criteria)
end

function bandpass!(gains,
                   data::Array{Complex64,3},
                   model::Array{Complex64,3},
                   interferometer::Interferometer,
                   criteria::StoppingCriteria)
    # Transpose the data and model arrays to create a better memory access pattern
    data  = permutedims(data, (3,1,2))
    model = permutedims(model,(3,1,2))

    workspace = Workspace_bandpass(interferometer)

    for β = 1:interferometer.Nfreq
        # Calibrate the X polarization
        bandpass_onechannel!(slice(gains,:,1,β),
                             slice( data,:,1,β),
                             slice(model,:,1,β),
                             interferometer,
                             workspace,
                             criteria)
        # Calibrate the Y polarization
        bandpass_onechannel!(slice(gains,:,2,β),
                             slice( data,:,4,β),
                             slice(model,:,4,β),
                             interferometer,
                             workspace,
                             criteria)
    end
    nothing
end

# Pre-allocate all the arrays needed by the calibration routine

immutable InnerWorkspace_bandpass
    data::Matrix{Complex64}
    model::Matrix{Complex64}
    normalization::Vector{Float64}
end

immutable Workspace_bandpass
    data::Matrix{Complex64}  # Hermitian matrix
    model::Matrix{Complex64} # Hermitian matrix
    gains::Vector{Complex64}
    oldgains::Vector{Complex64}

    rkworkspace::RKWorkspace{Complex64,4}
    innerworkspace::InnerWorkspace_bandpass
end

function Workspace_bandpass(interferometer::Interferometer)
    Nant = interferometer.Nant - length(interferometer.flaggedantennas)
    data  = Array(Complex64,Nant,Nant)
    model = Array(Complex64,Nant,Nant)
    gains = Array(Complex64,Nant)
    oldgains = Array(Complex64,Nant)
    normalization = Array(Float64,Nant)
    Workspace_bandpass(data,model,gains,oldgains,
                       RKWorkspace(gains,4),
                       InnerWorkspace_bandpass(data,model,normalization))
end

@doc """
Calibrate the complex gains from a single frequency channel using a two step process:

1. Get an initial estimate of the gains.
2. Iteratively improve that estimate.
""" ->
function bandpass_onechannel!(gains, data, model,
                              interferometer::Interferometer,
                              workspace::Workspace_bandpass,
                              criteria::StoppingCriteria)
    # Pack the visibilities into square, Hermitian matrices
    makesquare!(workspace.data,  data,interferometer.flaggedantennas)
    makesquare!(workspace.model,model,interferometer.flaggedantennas)

    # Take a first guess at the complex gains
    firstguess!(workspace)

    # Refine the estimate of the gains
    iter = 0
    converged = false
    while !converged && iter < criteria.maxiter
        workspace.oldgains[:] = workspace.gains
        rkstep!(workspace.gains,bandpass_step,workspace.innerworkspace,workspace.rkworkspace)
        if vecnorm(workspace.gains-workspace.oldgains)/vecnorm(workspace.oldgains) < criteria.tolerance
            converged = true
        end
        iter += 1
    end

    # Output
    idx = 1
    for ant = 1:interferometer.Nant
        ant in interferometer.flaggedantennas && continue
        gains[ant] = workspace.gains[idx]
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
function firstguess!(workspace::Workspace_bandpass)
    G = workspace.data./workspace.model
    λ,v = eigs(G,nev=1,which=:LM)
    w = sqrt(λ[1])
    @devec workspace.gains[:] = v.*w
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
function bandpass_step!(output,input,workspace)
    N = length(input)
    for i = 1:N
        output[i] = complex(0.)
        workspace.normalization[i] = 0.
    end
    for j = 1:N, i = 1:N
        i == j && continue
        z = conj(input[j])*workspace.model[i,j]
        output[i] += conj(z)*workspace.data[i,j]
        workspace.normalization[i] += abs2(z)
    end
    for i = 1:N
        output[i] = output[i]/workspace.normalization[i] - input[i]
    end
    nothing
end
const bandpass_step = RKInnerStep{:bandpass_step!}()

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

function χ2(workspace::Workspace_bandpass)
    χ2(workspace.data,workspace.model,workspace.gains)
end

