# Copyright (c) 2015 Michael Eastwood
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
    genvis(meta::Metadata, sources)

Generate model visibilities for the given list of sources.
"""
function genvis(meta::Metadata, sources)
    model = zeros(JonesMatrix, Nbase(meta), Nfreq(meta))
    flags = zeros(       Bool, Nbase(meta), Nfreq(meta))
    for source in sources
        genvis_onesource!(model, meta, source)
    end
    Visibilities(model, flags)
end

genvis(meta::Metadata, source::Source) = genvis(meta, [source])

function frame_and_phase_center(meta::Metadata)
    frame = reference_frame(meta)
    phase_center = measure(frame, meta.phase_center, dir"ITRF")
    frame, phase_center
end

"""
    refract_and_corrupt(meta, frame, source)

This function returns the flux and ITRF direction to the source
at each frequency channel after refraction through the ionosphere
and corruption by the primary beam.
"""
function refract_and_corrupt(meta, frame, source)
    pos  = position(meta)
    azel = measure(frame, source.direction, dir"AZEL")
    flux = zeros(HermitianJonesMatrix, Nfreq(meta))
    itrf = fill(Direction(dir"ITRF", 1.0, 0.0, 0.0), Nfreq(meta))
    for β = 1:Nfreq(meta)
        frequency = meta.channels[β]
        # refraction through the ionosphere
        # TODO: allow the plasma frequency to be set to something nonzero
        refracted_direction = refract(azel, frequency, 0.0)
        # corruption by the primary beam
        my_flux = (source.spectrum(frequency) |> linear) :: HermitianJonesMatrix
        jones   = meta.beam(frequency, refracted_direction) :: JonesMatrix
        flux[β] = congruence_transform(jones, my_flux)
        # convert the source direction into the correct coordinate system
        itrf[β] = azel_to_itrf(refracted_direction, pos)
    end
    flux, itrf
end

"""
    additional_precomputation(meta, frame, source)

This function is used to do any computations we would like to only do
once per source. For example it is used to compute the major and minor
axes for Gaussian sources.
"""
additional_precomputation(meta, frame, source) = nothing

"""
    baseline_coherency(source, frequency, antenna1, antenna2, variables)

This function computes the coherency between the two given antennas for
the given source. For a point source this is unity for all baselines.
For a Gaussian source this is a Gaussian function of the baseline length
and orientation.

Note that the `variables` argument is the output of the
`additional_precomputation` function.
"""
baseline_coherency(source, frequency, antenna1, antenna2, variables) = 1

function genvis_onesource!(visibilities, meta, source)
    frame, phase_center = frame_and_phase_center(meta)
    flux, itrf_direction = refract_and_corrupt(meta, frame, source)
    variables = additional_precomputation(meta, frame, source)
    for β = 1:Nfreq(meta)
        frequency = meta.channels[β]
        visibilities_onechannel = slice(visibilities, :, β)
        genvis_onesource_onechannel!(visibilities_onechannel, meta, source, frequency,
                                     flux[β], itrf_direction[β], phase_center, variables)
    end
    visibilities
end

function genvis_onesource!(visibilities, meta, source::MultiSource)
    for component in source.components
        genvis_onesource!(visibilities, meta, component)
    end
    visibilities
end

function genvis_onesource_onechannel!(visibilities, meta, source, frequency,
                                      flux, itrf, phase_center, variables)
    delays  = geometric_delays(meta.antennas, itrf, phase_center)
    fringes = delays_to_fringes(delays, frequency)
    for α = 1:Nbase(meta)
        idx1 = meta.baselines[α].antenna1
        idx2 = meta.baselines[α].antenna2
        antenna1 = meta.antennas[idx1]
        antenna2 = meta.antennas[idx2]
        fringe = fringes[idx1] * conj(fringes[idx2])
        coherency = baseline_coherency(source, frequency, antenna1, antenna2, variables)
        visibilities[α] += flux * fringe * coherency
    end
    visibilities
end

function additional_precomputation(meta, frame, source::GaussianSource)
    j2000 = measure(frame, source.direction, dir"J2000")
    rhat  = [j2000.x, j2000.y, j2000.z] # the source location
    north = [0, 0, 1] - j2000.z*rhat    # local north on the celestial sphere
    east  = cross(north, rhat)          # local east on the celestial sphere
    north = north / norm(north)
    east  = east / norm(east)
    θ = source.position_angle
    major_axis =  cos(θ)*north + sin(θ)*east
    minor_axis = -sin(θ)*north + cos(θ)*east
    major_width = π^2 * sin(source.major_fwhm)^2 / (4log(2))
    minor_width = π^2 * sin(source.minor_fwhm)^2 / (4log(2))
    # convert the major and minor axes to the ITRF coordinate system
    major_j2000 = Direction(dir"J2000", major_axis[1], major_axis[2], major_axis[3])
    minor_j2000 = Direction(dir"J2000", minor_axis[1], minor_axis[2], minor_axis[3])
    major_itrf = measure(frame, major_j2000, dir"ITRF")
    minor_itrf = measure(frame, minor_j2000, dir"ITRF")
    major_itrf, minor_itrf, major_width, minor_width
end

function baseline_coherency(source::GaussianSource, frequency, antenna1, antenna2, variables)
    λ = c / frequency
    u = (antenna1.position.x - antenna2.position.x) / λ
    v = (antenna1.position.y - antenna2.position.y) / λ
    w = (antenna1.position.z - antenna2.position.z) / λ
    major_axis, minor_axis, major_width, minor_width = variables
    # project the baseline onto the major and minor axes
    major_proj = u*major_axis.x + v*major_axis.y + w*major_axis.z
    minor_proj = u*minor_axis.x + v*minor_axis.y + w*minor_axis.z
    #if frequency == 24.768e6
    #end
    # flux is attenuated on long baselines due to souce being a Gaussian
    exp(-major_width*major_proj^2 - minor_width*minor_proj^2)
end

"""
    geometric_delays(antennas, source_direction::Direction, phase_center)

Compute the geometric delay to each antenna for a source in the far field
of the interferometer.
"""
function geometric_delays(antennas, source_direction::Direction, phase_center)
    # far field
    l = source_direction.x - phase_center.x
    m = source_direction.y - phase_center.y
    n = source_direction.z - phase_center.z
    delays = zeros(length(antennas))
    for i = 1:length(antennas)
        antenna_position = antennas[i].position
        x = antenna_position.x
        y = antenna_position.y
        z = antenna_position.z
        delays[i] = (x*l + y*m + z*n) / c
    end
    delays
end

"""
    geometric_delays(antennas, source_position::Position, phase_center)

Compute the geometric delay to each antenna for a source in the near field
of the interferometer.
"""
function geometric_delays(antennas, source_position::Position, phase_center)
    # near field
    l = -phase_center.x
    m = -phase_center.y
    n = -phase_center.z
    ξ = source_position.x
    η = source_position.y
    ζ = source_position.z
    D = sqrt(ξ^2 + η^2 + ζ^2)
    delays = zeros(length(antennas))
    for i = 1:length(antennas)
        antenna_position = antennas[i].position
        x = antenna_position.x
        y = antenna_position.y
        z = antenna_position.z
        delays[i] = (D - sqrt((x-ξ)^2 + (y-η)^2 + (z-ζ)^2) + x*l + y*m + z*n) / c
    end
    delays
end

doc"""
    delays_to_fringes(delays, frequency)

Compute
\[
    \exp\left(2\pi i \nu \tau\right)
\]
for each delay $\tau$ and the frequency $\nu$.
"""
function delays_to_fringes(delays, frequency)
    i2π = 1im * 2π
    Nant = length(delays)
    fringes = zeros(Complex128, Nant)
    for i = 1:Nant
        ϕ = i2π * frequency * delays[i]
        fringes[i] = exp(ϕ)
    end
    fringes
end

