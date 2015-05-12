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
    time   = Epoch(Measures.UTC,Quantity(ms["TIME",1],Second))
    antpos = slice(Table(ms[kw"ANTENNA"])["POSITION"],:,1)
    pos    = Measures.from_xyz_in_meters(Measures.ITRF,antpos[1],antpos[2],antpos[3])
    set!(frame,time)
    set!(frame,pos)
    frame
end

function ants(ms::Table)
    ant1 = ms["ANTENNA1"]
    ant2 = ms["ANTENNA2"]
    # (the +1 converts to a 1-based indexing scheme)
    ant1+1,ant2+1
end

function phase_dir(ms::Table)
    field = Table(ms[kw"FIELD"])
    dir = field["PHASE_DIR"]
    Direction(Measures.J2000,Quantity(dir[1],Radian),Quantity(dir[2],Radian))
end

@doc """
Get the flags from the dataset, but this information is stored in multiple locations.
Unify all these flags before returning.
""" ->
function flags(ms::Table)
    data_flags = ms["FLAG"]
    row_flags  = ms["FLAG_ROW"]
    for α = 1:length(row_flags)
        if row_flags[α]
            data_flags[:,:,α] = true
        end
    end
    data_flags
end

@doc """
Get the CORRECTED_DATA column if it exists. Otherwise settle for the DATA column.
""" ->
function corrected_data(ms::Table)
    if Tables.checkColumnExists(ms,"CORRECTED_DATA")
        return ms["CORRECTED_DATA"]
    else
        return ms["DATA"]
    end
end

