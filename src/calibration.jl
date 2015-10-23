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
    applycal!(ms::MeasurementSet, calibration::Calibration;
              apply_to_corrected = false, force_imaging_columns = false)

Apply the calibration to the given measurement set.

**Arguments:**

* `ms` - the measurement set to which the calibration will be applied
* `calibration` - the calibration that will be applied

**Keyword Arguments:**

* `apply_to_corrected` - if this is set to true, the calibration will be
    applied to the CORRECTED_DATA column instead of the DATA column
* `force_imaging_columns` - if this is set to true, the calibrated data
    will be written to the CORRECTED_DATA column regardless of whether
    or not the column already exists
"""
function applycal!(ms::MeasurementSet,
                   calibration::Calibration;
                   apply_to_corrected::Bool = false,
                   force_imaging_columns::Bool = false)
    data  = apply_to_corrected? get_corrected_data(ms) : get_data(ms)
    flags = get_flags(ms)
    applycal!(data,flags,calibration,ms.ant1,ms.ant2)
    set_corrected_data!(ms,data,force_imaging_columns)
    set_flags!(ms,flags)
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
    corrupt!(data::Array{Complex64,3}, cal::Calibration, ant1, ant2)

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

