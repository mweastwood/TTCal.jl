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

const c = 2.99792e+8 * Meter/Second

# These functions are necessary because type inference
# currently fails such that [1Meter, 2Meter] != [1, 2]*Meter

@doc """
A convenience function for attaching units to an array.
""" ->
function addunits{T}(array::Array{T},unit::SIUnits.SIUnit)
    output = Array(quantity(T,unit),size(array))
    for i = 1:length(array)
        output[i] = array[i]*unit
    end
    output
end

@doc """
A convenience function for stripping units from an array.
""" ->
function stripunits(array::Array)
    output = Array(typeof(array[1].val),size(array))
    for i = 1:length(array)
        output[i] = array[i].val
    end
    output
end

