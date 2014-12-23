@doc """
This type stores information about the interferometer that is not read from
a measurement set. It is intended to supplement the information readily
obtained from a measurement set.

# Fields
* `Nant`:
The number of antennas in the array.
* `Nfreq`:
The number of frequency channels within a sub-band.
* `refant`:
The reference antenna (whose phase will be fixed to zero). Because we do
not (yet) do any kind of polarization calibration, the x and y phases of
the reference antenna are both zeroed out.
* `flaggedantennas`:
A list of the antenna numbers that should be flagged by the calibration
routine. The first antenna is number 1.
* `delays`:
Contains the delay (in seconds) to each antenna (both polarizations).
This is used to estimate the decorrelation due to signals arriving at
the correlator at different times. This field is an array of size Nant
by Npol.
* `fftlength`:
The number of samples used in the calculation of the FFT. This number
is somewhat complicated by the presence of the polyphase filter bank,
which effectively lengthens the FFT.
* `samplerate`:
The sample rate of the ADCs (in Hz).
""" ->
type Interferometer
    Nant::Int
    Nfreq::Int
    refant::Int
    flaggedantennas::Vector{Int}
    delays::Array{Float64,2}
    fftlength::Int
    samplerate::Float64
end

@doc """
Constructor for the Interferometer type that initializes fields to sensible
defaults for the OVRO LWA.
""" ->
function LWA()
    Interferometer(256,109,1,zeros(Bool,256),zeros(Float64,256,2),8196,200e6)
end

