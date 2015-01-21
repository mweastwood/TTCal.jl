################################################################################
# Public Interface

@doc """
Fit the visibilities to a model of point sources. The input model needs
to have the positions of the point sources relatively close, but the flux
can be wildly off.
""" ->
function fitvis(ms::Table,
                sources::Vector{Source},
                criteria::StoppingCriteria)
    frame = ReferenceFrame()
    set!(frame,Epoch("UTC",ms["TIME",1]*Second))
    set!(frame,Measures.observatory(frame,"OVRO_MMA"))
    data = ms["CORRECTED_DATA"]
    flags = ms["FLAG"]
    uvw  = ms["UVW"]
    u = addunits(squeeze(uvw[1,:],1),Meter)
    v = addunits(squeeze(uvw[2,:],1),Meter)
    w = addunits(squeeze(uvw[3,:],1),Meter)

    spw = Table(ms[kw"SPECTRAL_WINDOW"])
    ν = addunits(spw["CHAN_FREQ",1],Hertz)

    fitvis(frame,data,flags,u,v,w,ν,sources,criteria)
end

function fitvis(frame::ReferenceFrame,
                data::Array{Complex64,3},
                flags::Array{Bool,3},
                u,v,w,ν,
                sources::Vector{Source},
                criteria::StoppingCriteria)
    output = abovehorizon(frame,sources)
    workspace = Workspace_fitvis(frame,data,flags,u,v,w,ν,output)
    outer_fitvis(workspace,criteria)
    output
end

################################################################################
# Workspace Definition

type Workspace_fitvis
    frame::ReferenceFrame
    data::Array{Complex64,3}
    flags_T::Array{Bool,3}
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

    fitmat::Matrix{Complex128}
    Δdata::Vector{Complex128}
end

function Workspace_fitvis(frame,data,flags,u,v,w,ν,sources)
    model = similar(data)
    data_minus_model = similar(data)
    data_minus_model_T = permutedims(similar(data),(3,2,1))

    Nparams = 2 + length(ν)
    params = Array(Float64,Nparams)
    oldparams = similar(params)

    rkworkspace = RKWorkspace(similar(params),4)

    Nbase = length(u)
    fitmat = Array(Complex128,2Nbase,3)
    Δdata  = Array(Complex128,2Nbase)

    Workspace_fitvis(frame,data,permutedims(flags,(3,2,1)),u,v,w,ν,sources,
                     model,data_minus_model,data_minus_model_T,
                     params,oldparams,
                     rkworkspace,
                     fitmat,Δdata)
end

################################################################################
# Outer Methods (outside rkstep!)

@doc """
1. Sort the sources in order of decreasing flux.
2. For each source, subtract all other sources.
3. Fit for the source.
""" ->
function outer_fitvis(workspace,criteria)
    u = workspace.u
    v = workspace.v
    w = workspace.w
    ν = workspace.ν
    flags = workspace.flags_T

    # Flag short baselines
    for α = 1:length(u)
        if sqrt(u[α]^2+v[α]^2+w[α]^2) < 70Meter
            flags[α,:,:] = true
        end
    end

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
        @time rkstep!(params,fitvis_step,workspace,rkworkspace)

        @show (params[1] - oldparams[1]) (params[2] - oldparams[2]) vecnorm(params-oldparams)/vecnorm(oldparams)
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
    @show iter

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
    flags = workspace.flags_T

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

    @time for β = 1:Nfreq
        l_arr[β],m_arr[β],flux[β] = inner_fitvis_onechannel(slice(data,:,β,1),
                                                            slice(flags,:,β,1),
                                                            slice(fringe,:,β),
                                                            flux[β],l,m,u,v,w,ν[β],
                                                            workspace.fitmat,
                                                            workspace.Δdata)
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

function inner_fitvis_onechannel(data,flags,fringe,flux,l,m,u,v,w,ν,fitmat,Δdata)
    n = sqrt(1-l^2-m^2)
    λ = c/ν

    Nbase = length(u)
    Nflag = sum(flags)
    N = Nbase-Nflag
    # Columns contain: [∂V/∂l ∂V/∂m ∂V/∂S]
    # Where V is the contribution to the visibility of a single point source,
    #       l and m specify the position of the point source, and
    #       S is the flux of the point source
    idx = 1
    @inbounds for α = 1:Nbase
        flags[α] && continue
        # (note the factor of 0.5 accounts for the fact that we're only
        # considering a single polarization)
        model = 0.5*flux*fringe[α]
        Δdata[idx] = data[α] - model
        fitmat[idx,1] = complex(0.,2.0π*(u[α].val-w[α].val*l/n)/λ.val) * model
        fitmat[idx,2] = complex(0.,2.0π*(v[α].val-w[α].val*m/n)/λ.val) * model
        fitmat[idx,3] = 0.5 * fringe[α]

        # Constrain the complex conjugates as well
        # (this forces the parameters to be real)
        Δdata[N+idx] = conj(Δdata[idx])
        fitmat[N+idx,1] = conj(fitmat[idx,1])
        fitmat[N+idx,2] = conj(fitmat[idx,2])
        fitmat[N+idx,3] = conj(fitmat[idx,3])

        idx += 1
    end
    Δparam = real(sub(fitmat,1:2N,:)\sub(Δdata,1:2N,:))
    l += Δparam[1]
    m += Δparam[2]
    flux += Δparam[3]
    l,m,flux
end

