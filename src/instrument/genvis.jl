# Copyright (c) 2015-2017 Michael Eastwood
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

#"""
#    genvis(metadata, beam, sources)
#
#*Description*
#
#Generate model visibilities for the given list of sources.
#
#If a beam model is provided, the model visibilities will be given as if
#they were observed with antennas that have the given primary beam.
#
#Note that bandwidth smearing and time smearing are currently not
#accounted for by this routine.
#
#*Arguments*
#
#* `metadata` - the metadata describing the interferometer
#* `beam` - the antenna primary beam (assumed to be constant if not given)
#* `sources` - the list of sources to include in the sky model
#"""
function genvis(metadata::Metadata, beam::AbstractBeam, sky::SkyModel; polarization=Full)
    dataset = Dataset(metadata, polarization=polarization)
    genvis!(dataset, beam, sky)
    dataset
end

function genvis(metadata::Metadata, beam::AbstractBeam, source::Source; polarization=Full)
    dataset = Dataset(metadata, polarization=polarization)
    genvis!(dataset, beam, source)
    dataset
end

function genvis!(dataset::Dataset, beam::AbstractBeam, sky::SkyModel)
    for source in sky.sources
        genvis!(dataset, beam, source)
    end
    dataset
end

function genvis!(dataset::Dataset, beam::AbstractBeam, source::Source)
    for shape in source.shapes
        genvis!(dataset, beam, shape)
    end
    dataset
end

function genvis!(dataset::Dataset, beam::AbstractBeam, shape::AbstractShape)
    metadata = dataset.metadata
    frame = ReferenceFrame(metadata)
    precomputation = additional_precomputation(frame, shape)
    for (idx, time) in enumerate(metadata.times)
        set!(frame, time)
        if isabovehorizon(frame, shape)
            phase_center = measure(frame, metadata.phase_centers[idx], dir"ITRF")
            direction    = measure(frame, shape.direction, dir"ITRF")
            for (jdx, frequency) in enumerate(metadata.frequencies)
                flux = observe_through_beam(shape, beam, frame, frequency)
                visibilities = dataset[jdx, idx]
                genvis_onesource_onechannel!(visibilities, shape, frequency,
                                             flux, metadata.positions, direction, phase_center,
                                             precomputation)
            end
        end
    end
    dataset
end

#"RFI sources need some special treatment."
#function genvis_onesource!(visibilities, meta, beam, source::RFISource)
#    frame = reference_frame(meta)
#    phase_center = measure(frame, meta.phase_center, dir"ITRF")
#    itrf_position = measure(frame, source.position, pos"ITRF")
#    for β = 1:Nfreq(meta)
#        frequency = meta.channels[β]
#        flux = source.spectrum(frequency) |> HermitianJonesMatrix
#        visibilities_onechannel = view(visibilities, :, β)
#        genvis_onesource_onechannel!(visibilities_onechannel, meta, source, frequency,
#                                     flux, itrf_position, phase_center, ())
#    end
#    visibilities
#end

function genvis_onesource_onechannel!(visibilities, shape, frequency,
                                      flux, positions, direction, phase_center,
                                      precomputation)
    delays  = geometric_delays(positions, direction, phase_center)
    fringes = delays_to_fringes(delays, frequency)
    for antenna1 = 1:Nant(visibilities), antenna2 = antenna1:Nant(visibilities)
        fringe = fringes[antenna1] * conj(fringes[antenna2])
        position1 = positions[antenna1]
        position2 = positions[antenna2]
        coherency = baseline_coherency(shape, frequency, position1, position2, precomputation)
        accumulate!(visibilities, antenna1, antenna2, flux*fringe*coherency)
    end
    visibilities
end

function observe_through_beam(shape::AbstractShape, beam::AbstractBeam,
                              frame::ReferenceFrame, frequency)
    direction = measure(frame, shape.direction, dir"AZEL")
    azimuth   = longitude(direction)
    elevation =  latitude(direction)
    jones = beam(frequency, azimuth, elevation)
    flux  = shape.spectrum(frequency) |> HermitianJonesMatrix
    flux  = congruence_transform(jones, flux)
    flux
end

function accumulate!(visibilities::Visibilities{Full, T}, antenna1, antenna2, value) where {T}
    visibilities[antenna1, antenna2] += value
end

function accumulate!(visibilities::Visibilities{Dual, T}, antenna1, antenna2, value) where {T}
    visibilities[antenna1, antenna2] += DiagonalJonesMatrix(value.xx, value.yy)
end

function accumulate!(visibilities::Visibilities{XX, T}, antenna1, antenna2, value) where {T}
    visibilities[antenna1, antenna2] += value.xx
end

function accumulate!(visibilities::Visibilities{YY, T}, antenna1, antenna2, value) where {T}
    visibilities[antenna1, antenna2] += value.yy
end

#"""
#    refract_and_corrupt(meta, beam, frame, source)
#
#This function returns the flux and ITRF direction to the source
#at each frequency channel after refraction through the ionosphere
#and corruption by the primary beam.
#
#**TODO** the plasma frequency is currently hard-coded to zero
#"""
#function refract_and_corrupt(meta, beam, frame, source)
#    pos  = position(meta)
#    azel = measure(frame, source.direction, dir"AZEL")
#    flux = zeros(HermitianJonesMatrix, Nfreq(meta))
#    itrf = fill(Direction(dir"ITRF", 1.0, 0.0, 0.0), Nfreq(meta))
#    for β = 1:Nfreq(meta)
#        frequency = meta.channels[β]
#        # refraction through the ionosphere
#        refracted_direction = refract(azel, frequency, 0.0)
#        # corruption by the primary beam
#        my_flux = source.spectrum(frequency) |> HermitianJonesMatrix
#        jones   = beam(frequency, refracted_direction) :: JonesMatrix
#        flux[β] = congruence_transform(jones, my_flux)
#        # convert the source direction into the correct coordinate system
#        itrf[β] = azel_to_itrf(refracted_direction, pos)
#    end
#    flux, itrf
#end

"""
    additional_precomputation(meta, frame, source)

This function is used to do any computations we would like to only do
once per source. For example it is used to compute the major and minor
axes for Gaussian sources.
"""
additional_precomputation(frame, source) = nothing

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

#function time_smearing(itrf, frequency, antenna1, antenna2)
#    ω = 2π / 86164.09054 # angular rotation speed of the Earth
#    cosδ = cos(latitude(itrf)) # cosine declination
#    τ = 13 # hard coded to LWA integration time
#
#    # get the east-west component of the baseline
#    u, v, w = get_uvw(frequency, antenna1, antenna2)
#    up = [antenna1.position.x, antenna1.position.y, antenna1.position.z]
#    north = [0, 0, 1]
#    east = cross(north, up)
#    east /= norm(east)
#    baseline_east = abs(u*east[1] + v*east[2] + w*east[3])
#
#    sinc(ω * cosδ * baseline_east * τ)
#end

#function bandwidth_smearing(itrf, antenna1, antenna2)
#    λ = c / 24e3 # hard coded to LWA channel width
#    u = (antenna1.position.x - antenna2.position.x) / λ
#    v = (antenna1.position.y - antenna2.position.y) / λ
#    w = (antenna1.position.z - antenna2.position.z) / λ
#    delay = u*itrf.x + v*itrf.y + w*itrf.z
#    sinc(delay)
#end

"""
    geometric_delays(positions, direction, phase_center)

Compute the geometric delay to each antenna for a source in the far field of the interferometer.
"""
function geometric_delays(positions::Vector{Position},
                          direction::Direction, phase_center::Direction)
    # far field
    delays = zeros(typeof(1.0*u"s"), length(positions))
    for i = 1:length(positions)
        path_difference = dot(positions[i], direction - phase_center)
        delays[i] = path_difference / u"c"
    end
    delays
end

#"""
#    geometric_delays(antennas, source_position::Position, phase_center)
#
#Compute the geometric delay to each antenna for a source in the near field
#of the interferometer.
#"""
#function geometric_delays(antennas, source_position::Position, phase_center)
#    # near field
#    l = -phase_center.x
#    m = -phase_center.y
#    n = -phase_center.z
#    ξ = source_position.x
#    η = source_position.y
#    ζ = source_position.z
#    D = sqrt(ξ^2 + η^2 + ζ^2)
#    delays = zeros(length(antennas))
#    for i = 1:length(antennas)
#        antenna_position = antennas[i].position
#        x = antenna_position.x
#        y = antenna_position.y
#        z = antenna_position.z
#        delays[i] = (D - sqrt((x-ξ)^2 + (y-η)^2 + (z-ζ)^2) + x*l + y*m + z*n) / c
#    end
#    delays
#end

doc"""
    delays_to_fringes(delays, frequency)

Compute $\exp(2πiντ)$ for each delay $τ$ and the frequency $ν$.
"""
function delays_to_fringes(delays, frequency)
    twoπ = float(2π)
    Nant = length(delays)
    fringes = zeros(Complex128, Nant)
    for i = 1:Nant
        fringes[i] = cis(twoπ * frequency * delays[i])
    end
    fringes
end

function local_north_east(frame, direction)
    j2000 = measure(frame, direction, dir"J2000")
    north = gram_schmidt(Direction(dir"J2000", 0, 0, 1), j2000)
    east  = cross(north, j2000)
    north, east
end

#function get_uvw(frequency, antenna1, antenna2)
#    λ = c / frequency
#    u = (antenna1.position.x - antenna2.position.x) / λ
#    v = (antenna1.position.y - antenna2.position.y) / λ
#    w = (antenna1.position.z - antenna2.position.z) / λ
#    u, v, w
#end

# TODO move functions below this comment out of this file

function additional_precomputation(frame, shape::Gaussian)
    north, east = local_north_east(frame, shape.direction)
    θ = shape.position_angle
    major_axis =  cos(θ)*north + sin(θ)*east
    minor_axis = -sin(θ)*north + cos(θ)*east
    major_width = π^2 * sin(shape.major_fwhm)^2 / (4log(2))
    minor_width = π^2 * sin(shape.minor_fwhm)^2 / (4log(2))
    # convert the major and minor axes to the ITRF coordinate system
    major_itrf = measure(frame, major_axis, dir"ITRF")
    minor_itrf = measure(frame, minor_axis, dir"ITRF")
    major_itrf, minor_itrf, major_width, minor_width
end

function baseline_coherency(::Gaussian, frequency, position1, position2, variables)
    major_axis, minor_axis, major_width, minor_width = variables
    # project the baseline onto the major and minor axes
    λ = u"c" / frequency
    baseline = position2 - position1
    major_proj = dot(baseline, major_axis) / λ
    minor_proj = dot(baseline, minor_axis) / λ
    # flux is attenuated on long baselines due to souce being a Gaussian
    exp(-major_width*major_proj^2 - minor_width*minor_proj^2)
end

#function baseline_coherency(source::DiskSource, frequency, antenna1, antenna2, variables)
#    u, v, w = get_uvw(frequency, antenna1, antenna2)
#    b = sqrt(u^2 + v^2 + w^2) # baseline length
#    δθ = source.radius
#    δθb = δθ*b
#    if δθb < eps(Float64)
#        return 1.0
#    else
#        return besselj1(2π*δθb)/(π*δθb)
#    end
#end
#
#function additional_precomputation(meta, frame, source::ShapeletSource)
#    north, east = local_north_east(frame, source.direction)
#    north_j2000 = Direction(dir"J2000", north[1], north[2], north[3])
#    east_j2000  = Direction(dir"J2000",  east[1],  east[2],  east[3])
#    north_itrf  = measure(frame, north_j2000, dir"ITRF")
#    east_itrf   = measure(frame,  east_j2000, dir"ITRF")
#    north_itrf, east_itrf
#end
#
#function baseline_coherency(source::ShapeletSource, frequency, antenna1, antenna2, variables)
#    u, v, w = get_uvw(frequency, antenna1, antenna2)
#    north, east = variables
#    β = source.scale
#    # project the baseline onto the sky
#    x = u*east.x  + v*east.y  + w*east.z
#    y = u*north.x + v*north.y + w*north.z
#    π2 = π*π
#    exponential = exp(-2π2*β^2*(x^2+y^2))
#    # compute the contribution from each shapelet component
#    out = zero(Complex128)
#    idx = 1
#    nmax = round(Int, sqrt(length(source.coeff))) - 1
#    for n2 = 0:nmax, n1 = 0:nmax
#        if source.coeff[idx] != 0
#            phase = (1im)^(n1+n2)
#            shapelet = phase * hermite(n1, 2π*x*β) * hermite(n2, 2π*y*β) * exponential
#            out += source.coeff[idx] * shapelet
#        end
#        idx += 1
#    end
#    out
#end

