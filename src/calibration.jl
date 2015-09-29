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

abstract Calibration

"""
Apply the calibration to the given measurement set.
"""
function applycal!(ms::Table,
                   calibration::Calibration;
                   apply_to_corrected::Bool = false,
                   force_imaging_columns::Bool = false)
    data  = apply_to_corrected? ms["CORRECTED_DATA"] : ms["DATA"]
    flags = MeasurementSets.flags(ms)
    ant1,ant2 = MeasurementSets.antennas(ms)
    applycal!(data,flags,calibration,ant1,ant2)
    if force_imaging_columns || Tables.checkColumnExists(ms,"CORRECTED_DATA")
        ms["CORRECTED_DATA"] = data
    else
        ms["DATA"] = data
    end
    ms["FLAG"] = flags
    data
end

function applycal!(data::Array{Complex64,3},
                   flags::Array{Bool,3},
                   cal::Calibration,
                   ant1,ant2)
    inverse_cal = invert(cal)
    corrupt!(data,flags,inverse_cal,ant1,ant2)
    data
end

"""
Corrupt the model data as if it had been observed
with an instrument with the given calibration.
"""
function corrupt!(data::Array{Complex64,3},
                  cal::Calibration,
                  ant1,ant2)
    flags = fill(false,size(data))
    corrupt!(data,flags,cal,ant1,ant2)
end

write(filename,calibration::Calibration) = JLD.save(filename,"cal",calibration)
read(filename) = JLD.load(filename,"cal")

include("ampcal.jl")
include("gaincal.jl")
include("polcal.jl")

