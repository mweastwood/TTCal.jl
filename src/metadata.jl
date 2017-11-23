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

struct Metadata
    frequencies   :: Vector{typeof{1.0*u"Hz"}}
    times         :: Vector{Epoch}
    positions     :: Vector{Metadata}  # ITRF
    phase_centers :: Vector{Direction} # ITRF (one per time)
end

Nfreq(metadata::Metadata) = length(metadata.frequencies)
Ntime(metadata::Metadata) = length(metadata.times)
Nant(metadata::Metadata)  = length(metadata.positions)

"Number of baselines (counting auto-correlations)"
function Nbase(metadata::Metadata)
    N = Nant(metadata)
    (N*(N+1)) ÷ 2
end

function merge!(lhs::Metadata, rhs::Metadata; axis=:frequency)
    if axis == :frequency
        append!(lhs.frequencies, rhs.frequencies)
    elseif axis == :time
        append!(lhs.times, rhs.times)
        append!(lhs.phase_centers, rhs.phase_centers)
    end
    lhs
end

"Read metadata from the given measurement set."
function Metadata(ms::Table)
    frequencies   = read_frequencies(ms)
    times         = read_times(ms)
    positions     = read_positions(ms)
    phase_centers = read_phase_centers(ms)
    Metadata(frequencies, times, positions, phase_centers)
end

"Read frequency channels from the `SPECTRAL_WINDOW` subtable."
function read_frequencies(ms::Table)
    spw_table   = ms[kw"SPECTRAL_WINDOW"]
    frequencies = spw_table["CHAN_FREQ", 1] .* u"Hz"
    Tables.close(spw_table)
    frequencies
end

"Read the times from the main table."
function read_times(ms::Table)
    time = ms["TIME", 1]
    [Epoch(epoch"UTC", time*u"s")]
end

"Read antenna positions from the `ANTENNA` subtable."
function read_positions(ms::Table)
    antenna_table = ms[kw"ANTENNA"]
    xyz = antenna_table["POSITION"]
    positions = Position[]
    for i = 1:size(xyz, 2)
        x = xyz[1, i]
        y = xyz[2, i]
        z = xyz[3, i]
        position = Position(pos"ITRF", x, y, z)
        push!(positions, position)
    end
    Tables.close(antenna_table)
    positions
end

"Read the phase directions from the `PHASE_DIR` subtable."
function read_phase_centers(ms::Table)
    field_table = ms[kw"FIELD"]
    dir = field_table["PHASE_DIR"]
    Tables.close(field_table)
    [Direction(dir"J2000", dir[1]*u"rad", dir[2]*u"rad")]
end

"Construct the reference frame of the interferometer."
function Tables.ReferenceFrame(metadata)
    frame = ReferenceFrame()
    set!(frame, metadata.times[1])
    set!(frame, position(metadata))
    frame
end

"Get the interferometer location by averaging the antenna positions."
position(metadata::Metadata) = mean(metadata.positions)

#"Read baselines from the main table."
#function read_baselines(ms::Table)
#    ant1 = ms["ANTENNA1"]
#    ant2 = ms["ANTENNA2"]
#    [Baseline(ant1[α]+1, ant2[α]+1) for α = 1:length(ant1)]
#end
#
#immutable UVW
#    u :: Vector{Float64}
#    v :: Vector{Float64}
#    w :: Vector{Float64}
#end
#
#function UVW(meta::Metadata)
#    u = zeros(Nbase(meta))
#    v = zeros(Nbase(meta))
#    w = zeros(Nbase(meta))
#    for α = 1:Nbase(meta)
#        antenna1 = meta.baselines[α].antenna1
#        antenna2 = meta.baselines[α].antenna2
#        r1 = meta.antennas[antenna1].position
#        r2 = meta.antennas[antenna2].position
#        u[α] = r1.x - r2.x
#        v[α] = r1.y - r2.y
#        w[α] = r1.z - r2.z
#    end
#    UVW(u, v, w)
#end

