################################################################################
# Public Interface

@doc """
Fit the visibilities to a model of point sources. The input model needs
to have the positions of the point sources relatively close, but the flux
can be wildly off.
""" ->
function fitvis(interferometer::Interferometer,
                ms::Table,
                sources::Vector{Source},
                criteria::StoppingCriteria)
    # TODO
end

function fitvis(frame::ReferenceFrame,
                data::Array{Complex64,3},
                u,v,w,ν,
                sources::Vector{Source},
                criteria::StoppingCriteria)
    output = [source for source in sources]
    workspace = Workspace_fitvis(frame,data,u,v,w,ν,output)
    outer_fitvis(workspace,criteria)
    output
end

################################################################################
# Workspace Definition

type Workspace_fitvis
    frame::ReferenceFrame
    data::Array{Complex64,3}
    u::Vector{quantity(Float64,Meter)}
    v::Vector{quantity(Float64,Meter)}
    w::Vector{quantity(Float64,Meter)}
    ν::Vector{quantity(Float64,Hertz)}
    sources::Vector{Source}

    model::Array{Complex64,3}
    data_minus_model::Array{Complex64,3}
    data_minus_model_T::Array{Complex64,3}

    params::Vector{Float64}
    oldparams::Vector{Float64}

    rkworkspace::RKWorkspace{Float64,4}
end

function Workspace_fitvis(frame,data,u,v,w,ν,sources)
    model = similar(data)
    data_minus_model = similar(data)
    data_minus_model_T = permutedims(similar(data),(3,2,1))

    Nparams = 2 + length(ν)
    params = Array(Float64,Nparams)
    oldparams = similar(params)

    rkworkspace = RKWorkspace(similar(params),4)

    Workspace_fitvis(frame,data,u,v,w,ν,sources,
                     model,data_minus_model,data_minus_model_T,
                     params,oldparams,
                     rkworkspace)
end

################################################################################
# Outer Methods (outside rkstep!)

@doc """
1. Sort the sources in order of decreasing flux.
2. For each source, subtract all other sources.
3. Fit for the source.
""" ->
function outer_fitvis(workspace,criteria)
    sources = workspace.sources
    data = workspace.data
    data_minus_model = workspace.data_minus_model
    data_minus_model_T = workspace.data_minus_model_T

    # 1. Sort the sources in order of decreasing flux.
    sort!(sources, by=source->source.flux, rev=true)
    # 2. For each source, subtract all other sources.
    for i = 1:length(sources)
        data_minus_model[:] = data
        for j = 1:length(sources)
            i == j && continue
            subtract_source!(workspace,sources[j])
        end
        data_minus_model_T[:] = permutedims(data_minus_model,(3,2,1))
        # 3. Fit for the source.
        sources[i] = outer_fitvis_singlesource(sources[i],workspace,criteria)
    end
end

function subtract_source!(workspace,source)
    frame = workspace.frame
    model = workspace.model
    data_minus_model = workspace.data_minus_model
    u = workspace.u
    v = workspace.v
    w = workspace.w
    ν = workspace.ν

    model[:] = visibilities(frame,source,u,v,w,ν)
    @devec data_minus_model[:] = data_minus_model-model
end

@doc """
1. Pack all the fitting paramaters into a vector.
2. Take a Runge-Kutta step until the parameters have converged.
3. Get the J2000 position of the source.
4. Fit a 2-component power law spectrum to the source.
""" ->
function outer_fitvis_singlesource(source,workspace,criteria)
    frame = workspace.frame
    ν = workspace.ν
    params = workspace.params
    oldparams = workspace.oldparams
    rkworkspace = workspace.rkworkspace

    # 1. Pack all the fitting paramaters into a vector.
    l,m  = getlm(frame,source)
    params[1] = l
    params[2] = m
    params[3:end] = getflux(source,ν)

    # 2. Take a Runge-Kutta step until the parameters have converged.
    iter = 0
    converged = false
    while !converged && iter < criteria.maxiter
        oldparams[:] = params
        rkstep!(params,fitvis_step,workspace,rkworkspace)

        # The first two elements are (l,m), the position of the source. These will
        # tend to be much smaller than the rest of the parameters, which is why I
        # test their convergence separately.
        # TODO: exclude (l,m) from the third check
        if (params[1] - oldparams[1]) < criteria.tolerance &&
           (params[2] - oldparams[2]) < criteria.tolerance &&
           vecnorm(params-oldparams)/vecnorm(oldparams) < criteria.tolerance
            converged = true
        end
        iter += 1
    end
    #error("stop")

    # 3. Get the J2000 position of the source.
    l = params[1]
    m = params[2]
    flux = sub(params,3:length(params))

    az = atan2(l,m)*Radian
    el = acos(sqrt(l.^2+m.^2))*Radian

    azel  = Direction("AZEL",az,el)
    j2000 = measure(frame,azel,"J2000")

    # 4. Fit a 2-component power law spectrum to the source.
    reffreq = 47e6*Hertz
    x = log10(stripunits(ν/reffreq))
    y = log10(flux)
    z = [ones(length(ν)) x x.^2]\y
    logflux   = z[1]
    index     = z[2:3]
    amplitude = 10.0.^logflux

    Source(source.name,j2000,amplitude,reffreq,index)
end

################################################################################
# Inner Methods (inside rkstep!)

@doc """
1. Generate the fringe pattern for all baselines, frequencies.
2. Fit for the position and flux of the source at each frequency.
3. Get a single position for the source.
""" ->
function inner_fitvis(output,input,workspace)
    u = workspace.u
    v = workspace.v
    w = workspace.w
    ν = workspace.ν
    data = workspace.data_minus_model_T

    l = input[1]
    m = input[2]
    flux = input[3:end]

    # 1. Generate the fringe pattern for all baselines, frequencies.
    fringe = permutedims(fringepattern(l,m,u,v,w,ν),(2,1))

    # fringe is ordered [baseline,frequency]
    #   data is ordered [baseline,frequency,polarization] (from an earlier transpose)

    # 2. Fit for the position and flux of the source at each frequency.
    Nfreq = length(ν)
    l_arr = Array(Float64,Nfreq)
    m_arr = Array(Float64,Nfreq)

    for β = 1:Nfreq
        l_arr[β],m_arr[β],flux[β] = inner_fitvis_onechannel(slice(data,:,β,1),
                                                            slice(fringe,:,β),
                                                            flux[β],l,m,u,v,w,ν[β])
    end

    # 3. Get a single position for the source.
    # (Weight the positions by ν^2 because resolution ~ 1/ν)
    weights = ν.^2
    l = sum(weights .* l_arr) ./ sum(weights)
    m = sum(weights .* m_arr) ./ sum(weights)

    output[1] = l-input[1]
    output[2] = m-input[2]
    output[3:end] = flux-input[3:end]
end
const fitvis_step = RKInnerStep{:inner_fitvis}()

function inner_fitvis_onechannel(data,fringe,flux,l,m,u,v,w,ν)
    n = sqrt(1-l^2-m^2)
    λ = c/ν

    Nbase  = length(u)
    fitmat = Array(Complex128,2Nbase,3)
    # Columns contain: [∂V/∂l ∂V/∂m ∂V/∂S]
    # Where V is the contribution to the visibility of a single point source,
    #       l and m specify the position of the point source, and
    #       S is the flux of the point source
    for α = 1:Nbase
        # (note the factor of 0.5 accounts for the fact that we're only
        # considering a single polarization)
        fitmat[α,1] = complex(0.,2.0π*(u[α]-w[α]*l/n)/λ) * 0.5flux * fringe[α]
        fitmat[α,2] = complex(0.,2.0π*(v[α]-w[α]*m/n)/λ) * 0.5flux * fringe[α]
        fitmat[α,3] = 0.5 * fringe[α]
        # Constrain the complex conjugates as well
        # (this forces the parameters to be real)
        fitmat[Nbase+α,1] = conj(fitmat[α,1])
        fitmat[Nbase+α,2] = conj(fitmat[α,2])
        fitmat[Nbase+α,3] = conj(fitmat[α,3])
    end
    model  = 0.5*flux*fringe
    Δdata  = vcat(data,conj(data)) - vcat(model,conj(model))
    Δparam = real(fitmat\Δdata)
    l += Δparam[1]
    m += Δparam[2]
    flux += Δparam[3]
    l,m,flux
end

