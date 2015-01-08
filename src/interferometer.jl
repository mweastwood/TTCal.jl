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
""" ->
type Interferometer
    Nant::Int
    Nfreq::Int
    refant::Int
    flaggedantennas::Vector{Int}
end

