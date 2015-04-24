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

# These are helper routines for interfacing with measurement sets

function uvw(ms::Table)
    uvw_arr = ms["UVW"]
    u = squeeze(uvw_arr[1,:],1)
    v = squeeze(uvw_arr[2,:],1)
    w = squeeze(uvw_arr[3,:],1)
    u,v,w
end

function freq(ms::Table)
    spw = Table(ms[kw"SPECTRAL_WINDOW"])
    ν = spw["CHAN_FREQ",1]
    ν
end

function reference_frame(ms::Table)
    frame = ReferenceFrame()
    set!(frame,Epoch(Measures.UTC,Quantity(ms["TIME",1],Second)))
    set!(frame,observatory("OVRO_MMA"))
    frame
end

function ants(ms::Table)
    ant1 = ms["ANTENNA1"]
    ant2 = ms["ANTENNA2"]
    # (the +1 converts to a 1-based indexing scheme)
    ant1+1,ant2+1
end

