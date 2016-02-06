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
    GenvisVariables

This type is simply a container for variables that need to be calculated
once and shared between functions while generating model visibilities.
"""
immutable GenvisVariables
    frame :: ReferenceFrame
    phase_center :: Direction # should always be ITRF
    function GenvisVariables(meta::Metadata)
        frame = reference_frame(meta)
        phase_center = measure(frame, meta.phase_center, dir"ITRF")
        new(frame, phase_center)
    end
end

"""
    genvis(meta::Metadata, sources)

Generate model visibilities for the given list of sources.
"""
function genvis(meta::Metadata, sources)
    variables = GenvisVariables(meta)
    model = zeros(JonesMatrix, Nbase(meta), Nfreq(meta))
    flags = zeros(       Bool, Nbase(meta), Nfreq(meta))
    for source in sources
        genvis_onesource!(model, meta, variables, source)
    end
    Visibilities(model, flags)
end

function genvis_onesource!(visibilities, meta, variables, source::PointSource)
    direction = measure(variables.frame, source.direction, dir"AZEL")
    phase_center = variables.phase_center
    for β = 1:Nfreq(meta)
        frequency = meta.channels[β]
        # refraction through the ionosphere
        # TODO: allow the plasma frequency to be set to something nonzero
        refracted_direction = refract(direction, frequency, 0.0)
        # corruption by the primary beam
        flux  = (source.spectrum(frequency) |> linear) :: HermitianJonesMatrix
        jones = meta.beam(frequency, refracted_direction) :: JonesMatrix
        flux  = congruence_transform(jones, flux)
        # convert the source direction into the correct coordinate system
        itrf_direction = azel_to_itrf(refracted_direction, position(meta))
        # generate the visibilities!
        visibilities_onechannel = slice(visibilities, :, β)
        genvis_onesource_onechannel!(visibilities_onechannel, meta,
                                     flux, itrf_direction, phase_center, frequency)
    end
    visibilities
end

function genvis_onesource!(visibilities, meta, variables, source::MultiSource)
    for component in source.components
        genvis_onesource!(visibilities, meta, variables, component)
    end
    visibilities
end

function genvis_onesource!(visibilities, meta, variables, source::RFISource)
    position = measure(variables.frame, source.position, pos"ITRF")
    phase_center = variables.phase_center
    for β = 1:Nfreq(meta)
        frequency = meta.channels[β]
        # we don't corrupt the flux with the primary beam here because
        # the spectrum of the RFI was likely measured through the primary beam
        flux  = (source.spectrum(frequency) |> linear) :: HermitianJonesMatrix
        # generate the visibilities!
        visibilities_onechannel = slice(visibilities, :, β)
        genvis_onesource_onechannel!(visibilities_onechannel, meta,
                                     flux, position, phase_center, frequency)
    end
    visibilities
end

function genvis_onesource_onechannel!(visibilities, meta, flux, source, phase_center, frequency)
    delays  = geometric_delays(meta.antennas, source, phase_center)
    fringes = delays_to_fringes(delays, frequency)
    for α = 1:Nbase(meta)
        antenna1 = meta.baselines[α].antenna1
        antenna2 = meta.baselines[α].antenna2
        fringe = fringes[antenna1] * conj(fringes[antenna2])
        visibilities[α] += flux * fringe
    end
    visibilities
end

function geometric_delays(antennas, source_direction::Direction, phase_center)
    # farfield
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

function geometric_delays(antennas, source_position::Position, phase_center)
    # nearfield
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

