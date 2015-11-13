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
    immutable MeasurementSet

This type is a wrapper around `CasaCore.Tables.Table` that
is intended to simplify most of the common interactions
between TTCal and measurement sets.
"""
immutable MeasurementSet
    table::Table
    frame::ReferenceFrame
    phase_direction::Direction{dir"J2000"}
    u::Vector{Float64}
    v::Vector{Float64}
    w::Vector{Float64}
    ν::Vector{Float64}
    ant1::Vector{Int32}
    ant2::Vector{Int32}
    Nant::Int
    Nbase::Int
    Nfreq::Int
end

"""
    MeasurementSet(name)

Open the measurement set at the given location. Assorted
quantities that are commonly used by TTCal are automatically
loaded and stored in fields.
"""
function MeasurementSet(name)
    ms = Table(name)
    antenna_table = Table(ms[kw"ANTENNA"])
    field_table   = Table(ms[kw"FIELD"])
    spw_table     = Table(ms[kw"SPECTRAL_WINDOW"])

    frame = ReferenceFrame()
    time  = ms["TIME",1]
    position = antenna_table["POSITION",1]
    set!(frame,Epoch(epoch"UTC",Quantity(time,"s")))
    set!(frame,Position(pos"ITRF",position[1],position[2],position[3]))

    dir = field_table["PHASE_DIR"]
    phase_direction = Direction(dir"J2000",Quantity(dir[1],"rad"),
                                           Quantity(dir[2],"rad"))

    uvw_arr = ms["UVW"]
    u = squeeze(uvw_arr[1,:],1)
    v = squeeze(uvw_arr[2,:],1)
    w = squeeze(uvw_arr[3,:],1)

    ν = spw_table["CHAN_FREQ",1]

    # the +1 converts to a 1-based indexing scheme
    ant1 = ms["ANTENNA1"] + 1
    ant2 = ms["ANTENNA2"] + 1

    Nant = numrows(antenna_table)
    Nbase = length(u)
    Nfreq = length(ν)

    unlock(antenna_table)
    unlock(field_table)
    unlock(spw_table)

    MeasurementSet(ms,frame,phase_direction,
                   u,v,w,ν,ant1,ant2,
                   Nant,Nbase,Nfreq)
end

Tables.unlock(ms::MeasurementSet) = Tables.unlock(ms.table)

"""
    get_flags(ms::MeasurementSet)

Get the flags from the dataset, but this information is stored in multiple locations.
Unify all these flags before returning.
"""
function get_flags(ms::MeasurementSet)
    data_flags = ms.table["FLAG"]
    row_flags  = ms.table["FLAG_ROW"]
    for α = 1:length(row_flags)
        if row_flags[α]
            data_flags[:,:,α] = true
        end
    end
    data_flags
end

function set_flags!(ms::MeasurementSet,flags)
    ms.table["FLAG"] = flags
end

function get_data(ms::MeasurementSet)
    ms.table["DATA"]
end

function set_model_data!(ms::MeasurementSet, model,
                         force_imaging_columns = false)
    if force_imaging_columns || Tables.checkColumnExists(ms.table,"MODEL_DATA")
        ms.table["MODEL_DATA"] = model
    end
end

"""
    get_corrected_data(ms::MeasurementSet)

Get the CORRECTED_DATA column if it exists. Otherwise settle for the DATA column.
"""
function get_corrected_data(ms::MeasurementSet)
    if Tables.checkColumnExists(ms.table,"CORRECTED_DATA")
        return ms.table["CORRECTED_DATA"]
    else
        return ms.table["DATA"]
    end
end

function set_corrected_data!(ms::MeasurementSet, data,
                             force_imaging_columns = false)
    if force_imaging_columns || Tables.checkColumnExists(ms.table,"CORRECTED_DATA")
        ms.table["CORRECTED_DATA"] = data
    else
        ms.table["DATA"] = data
    end
end

"""
    flag_short_baselines!(flags, minuvw, u, v, w, ν)

Flag all of the baselines whose length is less than `minuvw` wavelengths.

This is a common operation that can mitigate contamination by unmodeled
diffuse emission.
"""
function flag_short_baselines!(flags, minuvw, u, v, w, ν)
    Nbase = length(u)
    Nfreq = length(ν)
    for α = 1:Nbase, β = 1:Nfreq
        if u[α]^2 + v[α]^2 + w[α]^2 < (minuvw*c/ν[β])^2
            flags[:,β,α] = true
        end
    end
    flags
end

