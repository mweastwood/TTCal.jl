@doc """
Fit the visibilities to a model of point sources. The input model needs
to have the positions of the point sources relatively close, but the flux
can be wildly off. This function returns a list of the fitted fluxes
and positions at each frequency.
""" ->
function fitvisibilities(frame::ReferenceFrame,
                         data::Array{Complex64,3},
                         source::Source,
                         u,v,w,ν)
    Nfreq = length(ν)
    l,m = getlm(frame,source)
    fringe = fringepattern(l,m,u,v,w,ν)

    # Transpose the data and model arrays to improve
    # the memory access pattern.
    data   = permutedims(data, (3,2,1)) # now (baseline,frequency,polarization)
    fringe = permutedims(fringe, (2,1)) # now (baseline,frequency)

    # TODO: apply flags and baseline filtering

    #for β = 1:Nfreq
    for β = 1:1
        flux = getflux(source,ν[β])
        fitvisibilities(slice(data,:,β,1),
                        slice(fringe,:,β),
                        flux,l,m,u,v,w,ν[β])
    end
end

function fitvisibilities(data,fringe,flux,l,m,u,v,w,ν)
    n = sqrt(1-l^2-m^2)
    λ = c/ν

    Nbase  = length(u)
    fitmat = Array(Complex64,2Nbase,3)
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

