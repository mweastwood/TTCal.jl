################################################################################
# Public Interface

@doc """
This function extracts the spectrum in a given direction by means of an
inverse discrete Fourier transform. Note that no gridding is performed,
so this does *not* use a fast Fourier transform. However, the inverse
discrete Fourier transform *is* the least squares estimator for the flux
in a given direction (if all baselines are weighted equally).
""" ->
function getspec(ms::Table,
                 dir::Direction)
    data = ms["CORRECTED_DATA"]
    frame = reference_frame(ms)
    l,m = getlm(frame,dir)
    u,v,w = uvw(ms)
    ν = freq(ms)
    ant1,ant2 = ants(ms)
    flags = ms["FLAG"]
    getspec(data,l,m,u,v,w,ν,ant1,ant2,flags)
end

function getspec{T<:Complex}(data::Array{T,3},
                             l::Float64,
                             m::Float64,
                             u::Vector{quantity(Float64,Meter)},
                             v::Vector{quantity(Float64,Meter)},
                             w::Vector{quantity(Float64,Meter)},
                             ν::Vector{quantity(Float64,Hertz)},
                             ant1::Vector{Int32},
                             ant2::Vector{Int32},
                             flags::Array{Bool,3})
    fringe = fringepattern(l,m,u,v,w,ν)
    spec   = zeros(Float64,4,length(ν))
    count  = zeros(Int,1,length(ν)) # The number of baselines used in the calculation
    for α = 1:length(u)
        # Don't use auto-correlations
        ant1[α] == ant2[α] && continue
        # Don't use short baselines
        #sqrt(u[α]^2+v[α]^2+w[α]^2) > 10Meter && continue
        for β = 1:length(ν)
            any(flags[:,β,α]) && continue
            # Taking the real part of A is equivalent to
            # computing (A + conj(A))/2. The conjugate of A,
            # in this case, is the baseline going in the
            # opposite direction. Including this information
            # constrains the spectrum to be real.
            z = conj(fringe[β,α])
            spec[1,β] += real(data[1,β,α]*z)
            spec[2,β] += real(data[2,β,α]*z)
            spec[3,β] += real(data[3,β,α]*z)
            spec[4,β] += real(data[4,β,α]*z)
            count[β] += 1
        end
    end
    spec = spec./count
    spec
end

