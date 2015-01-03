@doc """
Fit the visibilities to a model of point sources. The input model needs
to have the positions of the point sources relatively close, but the flux
can be wildly off.
""" ->
function fitvisibilities(frame::ReferenceFrame,
                         data::Array{Complex64,3},
                         sources::Vector{Source},
                         u,v,w,ν)
    # TODO: Sort sources by flux before doing this.
    # The idea is that small errors on bright sources (looking at you Cyg A)
    # will propagate as irrecoverable large errors on fainter sources.
    # Therefore we want to make sure we get the bright sources right first.
    for i = 1:2
        sources_omitting_current = Source[]
        for j = 1:length(sources)
            i == j && continue
            push!(sources_omitting_current,sources[j])
        end
        model = visibilities(frame,sources_omitting_current,u,v,w,ν)
        data_minus_model = data-model
        newsource = fitvisibilities(frame,data_minus_model,
                                    sources[i],u,v,w,ν)
        sources[i] = newsource
    end
end

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

    lfit = Array(Float64,Nfreq)
    mfit = Array(Float64,Nfreq)
    fluxfit = Array(Float64,Nfreq)
    for β = 1:Nfreq
        flux = getflux(source,ν[β])
        l_,m_,flux_ = fitvisibilities(slice(data,:,β,1),
                                      slice(fringe,:,β),
                                      flux,l,m,u,v,w,ν[β])
        lfit[β] = l_
        mfit[β] = m_
        fluxfit[β] = flux_
    end
    fitsource(frame,source.name,ν,lfit,mfit,fluxfit)
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

function fitsource(frame,name,ν,l,m,flux)
    # Get the average position of each source
    # TODO: weight higher frequencies more (higher resolution)
    lnew = mean(l)
    mnew = mean(m)

    az = atan2(lnew,mnew)
    el = acos(sqrt(lnew.^2+mnew.^2))

    azel  = Direction("AZEL",Quantity(az,"rad"),Quantity(el,"rad"))
    j2000 = measure(frame,azel,"J2000")

    # Fit a 2-component spectrum to the source
    reffreq = 47e6
    x = log10(ν/reffreq)
    y = log10(flux)
    z = [ones(length(ν)) x x.^2]\y
    logflux   = z[1]
    index     = z[2:3]
    amplitude = 10.0.^logflux

    Source(name,j2000.m[1],j2000.m[2],amplitude,reffreq,index)
end

