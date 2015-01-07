@doc """
Fit the visibilities to a model of point sources. The input model needs
to have the positions of the point sources relatively close, but the flux
can be wildly off.
""" ->
function fitvis(frame::ReferenceFrame,
                data::Array{Complex64,3},
                sources::Vector{Source},
                u,v,w,ν)
    # TODO: Sort sources by flux before doing this.
    # The idea is that small errors on bright sources (looking at you Cyg A)
    # will propagate as irrecoverable large errors on fainter sources.
    # Therefore we want to make sure we get the bright sources right first.
    for i = 1:length(sources)
        sources_omitting_current = Source[]
        for j = 1:length(sources)
            i == j && continue
            push!(sources_omitting_current,sources[j])
        end
        model = visibilities(frame,sources_omitting_current,u,v,w,ν)
        data_minus_model = data-model
        sources[i] = fitvis(frame,data_minus_model,sources[i],u,v,w,ν)
    end
end

function fitvis(frame::ReferenceFrame,
                data::Array{Complex64,3},
                source::Source,
                u,v,w,ν)
    args = fitvis_args(permutedims(data,(3,2,1)),u,v,w,ν)
    params = [getlm(frame,source)..., getflux(source,ν)]

    x  = params
    x′ = similar(x)
    k  = [similar(x) for i = 1:4]

    maxiter = 30
    tol  = 1e-5
    iter = 0
    converged = false
    while !converged && iter < maxiter
        xold = copy(x)
        rkstep!(x,x′,k,args,fitvis_step,RK4)

        # TODO: implement a better stopping criterion that
        # considers the change in position (l,m) and change
        # in flux separately (because these typically have
        # two different orders of magnitude)
        if vecnorm(x-xold)/vecnorm(xold) < tol
            converged = true
        end
        iter += 1
    end

    l = x[1]
    m = x[2]
    flux = x[3:end]

    az = atan2(l,m)*Radian
    el = acos(sqrt(l.^2+m.^2))*Radian

    azel  = Direction("AZEL",az,el)
    j2000 = measure(frame,azel,"J2000")

    # Fit a 2-component spectrum to the source
    reffreq = 47e6*Hertz
    x = log10(stripunits(ν/reffreq))
    y = log10(flux)
    z = [ones(length(ν)) x x.^2]\y
    logflux   = z[1]
    index     = z[2:3]
    amplitude = 10.0.^logflux

    Source(source.name,j2000,amplitude,reffreq,index)
end

type fitvis_args
    data::Array{Complex64,3}
    u::Vector{quantity(Float64,Meter)}
    v::Vector{quantity(Float64,Meter)}
    w::Vector{quantity(Float64,Meter)}
    ν::Vector{quantity(Float64,Hertz)}
end

function fitvis_step!(output,input,args)
    l = input[1]
    m = input[2]
    flux = input[3:end]

    u = args.u
    v = args.v
    w = args.w
    ν = args.ν
    data = args.data
    fringe = permutedims(fringepattern(l,m,u,v,w,ν),(2,1))

    # fringe is ordered [baseline,frequency]
    #   data is ordered [baseline,frequency,polarization] (from an earlier transpose)

    Nfreq = length(ν)
    l_arr = Array(Float64,Nfreq)
    m_arr = Array(Float64,Nfreq)

    for β = 1:Nfreq
        l_arr[β],m_arr[β],flux[β] = fitvis_onechannel(slice(data,:,β,1),
                                                      slice(fringe,:,β),
                                                      flux[β],l,m,u,v,w,ν[β])
    end

    # Weight the positions by ν^2 because resolution ~ 1/ν
    weights = ν.^2
    l = sum(weights .* l_arr) ./ sum(weights)
    m = sum(weights .* m_arr) ./ sum(weights)

    output[1] = l-input[1]
    output[2] = m-input[2]
    output[3:end] = flux-input[3:end]
end
const fitvis_step = RKInnerStep{:fitvis_step!}()

function fitvis_onechannel(data,fringe,flux,l,m,u,v,w,ν)
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

