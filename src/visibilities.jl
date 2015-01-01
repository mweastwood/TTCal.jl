const c = 2.99792e+8 # m/s

# TODO
# - incorporate polarization
# - re-include the delay decorrelation

@doc """
Calculate model visibilities for placement in the MODEL_DATA column of the
given measurement set. No gridding is performed, so the runtime of this
naive algorithm scales os O(Nbase×Nsource)

# Inputs
* `interferometer`:
An instance of the `Interferometer` type that contains information
about the interferometer.
* `ms`:
A measurement set (opened as a CasaCore table).
* `sources`:
A list of `Source`s to use in the model.

# Outputs
* `model`:
The MODEL_DATA column of the measurement set.
""" ->
function visibilities(ms::Table,
                      sources::Vector{Source})
    uvw = ms["UVW"]
    u = uvw[1,:]
    v = uvw[2,:]
    w = uvw[3,:]
    spw = Table(ms[kw"SPECTRAL_WINDOW"])
    ν = spw["CHAN_FREQ",1]

    frame = ReferenceFrame()
    set!(frame,Epoch("UTC",Quantity(ms["TIME",1],"s")))
    set!(frame,Measures.observatory(frame,"OVRO_MMA"))
    
    visibilities(frame,sources,u,v,w,ν)
end

function visibilities(frame::ReferenceFrame,sources::Vector{Source},
                      u,v,w,ν)
    Nbase   = length(u)
    Nfreq   = length(ν)
    model = zeros(Complex64,4,Nfreq,Nbase)
    visibilities!(model,frame,sources,u,v,w,ν)
    model
end

function visibilities!(model::Array{Complex64,3},
                       frame::ReferenceFrame,sources::Vector{Source},
                       u,v,w,ν)
    Δν = ν[2] - ν[1]
    Nbase  = length(u)
    Nfreq  = length(ν)
    fringe = Array(Complex64,Nfreq)
    for source in sources
        # Get the position and flux of the source
        l,m = getlm(frame,source.ra,source.dec)
        n = sqrt(1-l^2-m^2)
        flux = getflux(source,ν)

        for α = 1:Nbase
            # Get the fringe pattern for the baseline
            τ = 2π*(u[α]*l+v[α]*m+w[α]*n)/c
            ϕ = τ*ν[1]
            Δϕ = τ*Δν
            fringepattern!(fringe,ϕ,Δϕ)

            # Calculate the contribution to the visibility
            for β = 1:Nfreq
                model[1,β,α] += 0.5*flux[β]*fringe[β] # xx
                model[4,β,α] += 0.5*flux[β]*fringe[β] # yy
            end
        end
    end
end

@doc """
Compute exp(i(ϕ+nΔϕ)) where ϕ and Δϕ define an equally space
grid of points where n = 1 to N.

Using the sine and cosine angle addition rules, you can define
an iterative method such that you only need to compute sines
and cosines for a single iteration.
""" ->
function fringepattern!{T<:Complex}(output::Array{T,1},ϕ,Δϕ)
    N = length(output)
    sin_Δϕ = sin(Δϕ)
    cos_Δϕ = cos(Δϕ)
    output[1] = complex(cos(ϕ),sin(ϕ))
    for n = 1:N-1
        output[n+1] = complex(real(output[n])*cos_Δϕ - imag(output[n])*sin_Δϕ,
                              imag(output[n])*cos_Δϕ + real(output[n])*sin_Δϕ)
    end
    nothing
end

function fringepattern(ϕ,Δϕ,N)
    output = Array(Complex128,N)
    fringepattern!(output,ϕ,Δϕ)
    output
end

