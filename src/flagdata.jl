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

################################################################################
# Public Interface

@doc """
Clear all of the flags in the measurement set.
""" ->
function clearflags!(ms::Table)
    N = numrows(ms)
    Nfreq = length(freq(ms))
    flags = zeros(Bool,4,Nfreq,N)
    ms["FLAG"] = flags
    flags
end

@doc """
Apply antenna flags to the measurement set.
""" ->
function flagdata!(ms::Table,flaggedantennas::Vector{Int})
    N = numrows(ms)
    Nfreq = length(freq(ms))
    ant1,ant2 = ants(ms)
    flags = ms["FLAG"]
    for α = 1:N
        if ant1[α] in flaggedantennas || ant2[α] in flaggedantennas
            flags[:,:,α] = true
        end
    end
    ms["FLAG"] = flags
    flags
end

