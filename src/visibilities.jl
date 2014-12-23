# TODO: Make visibility generation a lot faster!
const c = 2.99792e+8 # m/s

#=
@doc """
Calculates the visibility for a given point source on a given baseline.
This function accounts for the decorrelation due to input signals
arriving at the correlator at different times.

# Inputs
* `interferometer`:
An instance of the `Interferometer` type that contains information
about the interferometer.
* `ant1`, `ant2`:
The index of the antennas composing the given baseline. These are
used for accessing properties from within `interferometer`.
* `ν`:
The frequency in Hz.
* `u`, `v`, `w`:
The baseline (in the usual coordinate system) measured in units of
centimeters (not wavelengths).
* `l`, `m`:
The position of the source (in the usual coordinate system).
* `flux`:
The flux of the source. The source is assumed to be unpolarized.

# Outputs
* `vis_xx`: The xx visibility.
* `vis_yy`: The yy visibility.
""" ->
=#
#=
function visibility(interferometer::Interferometer,ant1,ant2,ν,u,v,w,l,m,flux)
    n = sqrt(1-l^2-m^2)
    geometric_delay = (u*l + v*m + w*n)/c
    signal_delay_x  = interferometer.delays[ant2,1] - interferometer.delays[ant1,1]
    signal_delay_y  = interferometer.delays[ant2,2] - interferometer.delays[ant1,2]
    decorrelation_x = 1-(abs(geometric_delay+signal_delay_x)
                        / (interferometer.fftlength/interferometer.samplerate))
    decorrelation_y = 1-(abs(geometric_delay+signal_delay_y)
                        / (interferometer.fftlength/interferometer.samplerate))

    vis = flux*exp(complex(0.,2π*ν*geometric_delay))
    0.5*vis*decorrelation_x,0.5*vis*decorrelation_y
end
=#

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
An instance of the `MeasurementSet` type.
* `sources`:
A list of `Source`s to use in the model.

# Outputs
* `model`:
The MODEL_DATA column of the measurement set.
""" ->
function visibilities(interferometer::Interferometer,
                      ms::MeasurementSet,
                      sources::Vector{Source})
    u,v,w = getUVW(ms)
    ν = getFreq(ms)
    Δν = ν[2] - ν[1]
    Nbase   = length(u)
    Nfreq   = length(ν)
    Nsource = length(sources)

    frame = ReferenceFrame()
    set!(frame,getTime(ms))
    set!(frame,observatory(frame,"OVRO_MMA"))

    fringe = Array(Complex64,Nfreq)
    model = zeros(Complex64,4,Nfreq,Nbase)
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
    model
end

@doc """
Compute exp(i(ϕ+nΔϕ)) where ϕ and Δϕ define an equally space
grid of points where n = 1 to N.

Using the sine and cosine angle addition rules, you can define
an iterative method such that you only need to compute sines
and cosines for a single iteration.
""" ->
function fringepattern!(output,ϕ,Δϕ)
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

# TODO: move this to sourcemodel.jl and rework
@doc """
Convert a given RA and dec to the standard radio coordinate system.
""" ->
function getlm(frame,ra,dec)
    dir  = Direction("J2000",ra,dec)
    azel = measure(frame,"AZEL",dir)
    az = azel.m[1].value
    el = azel.m[2].value
    l = cos(el)*sin(az)
    m = cos(el)*cos(az)
    l::Float64,m::Float64
end

