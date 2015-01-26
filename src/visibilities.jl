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
    u,v,w = uvw(ms)
    ν = freq(ms)
    frame = reference_frame(ms)
    visibilities(frame,sources,u,v,w,ν)
end

function visibilities(frame::ReferenceFrame,
                      sources::Vector{Source},
                      u,v,w,ν)
    model = zeros(Complex64,4,length(ν),length(u))
    visibilities!(model,frame,sources,u,v,w,ν)
    model
end

function visibilities!(model::Array{Complex64,3},
                       frame::ReferenceFrame,
                       sources::Vector{Source},
                       u,v,w,ν)
    for source in sources
        visibilities!(model,frame,source,u,v,w,ν)
    end
end

function visibilities(frame::ReferenceFrame,
                      source::Source,
                      u,v,w,ν)
    model = zeros(Complex64,4,length(ν),length(u))
    visibilities!(model,frame,source,u,v,w,ν)
    model
end

function visibilities!(model::Array{Complex64,3},
                       frame::ReferenceFrame,
                       source::Source,
                       u,v,w,ν)
    Δν = ν[2] - ν[1]
    Nbase = length(u)
    Nfreq = length(ν)

    # Get the position and flux of the source
    l,m = getlm(frame,source)
    n = sqrt(1-l^2-m^2)
    flux = getflux(source,ν)

    fringe = fringepattern(l,m,u,v,w,ν)
    for α = 1:Nbase, β = 1:Nfreq
        model[1,β,α] += 0.5*flux[β]*fringe[β,α] # xx
        model[4,β,α] += 0.5*flux[β]*fringe[β,α] # yy
    end
    nothing
end

function fringepattern(l,m,u,v,w,ν)
    output = Array(Complex64,length(ν),length(u))
    fringepattern!(output,l,m,u,v,w,ν)
    output
end

function fringepattern!{T<:Complex}(output::Array{T,2},l,m,u,v,w,ν)
    n = sqrt(1-l^2-m^2)
    Δν = ν[2] - ν[1]
    Nbase = length(u)
    Nfreq = length(ν)

    fringe = Array(Complex64,Nfreq)
    for α = 1:Nbase
        # Get the fringe pattern for the baseline
        τ = 2π*(u[α]*l+v[α]*m+w[α]*n)/c
        ϕ = τ*ν[1]
        Δϕ = τ*Δν
        fringepattern!(fringe,ϕ,Δϕ)

        # Calculate the contribution to the visibility
        for β = 1:Nfreq
            output[β,α] = fringe[β]
        end
    end
    nothing
end

function fringepattern(ϕ,Δϕ,N::Int)
    output = Array(Complex64,N)
    fringepattern!(output,ϕ,Δϕ)
    output
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

