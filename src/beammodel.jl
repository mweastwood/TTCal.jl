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
Generate beam model to attenuate source model flux (to be used in genvis.jl).
As an initial attempt, uses the naive (1-cos(el)^1.6) beam, normalized to 1 at zenith.
"""
function beammodel(frame,phase_dir,l,m)
    dir   = lm2dir(phase_dir,l,m)
    az,el = dir2azel(frame,dir)

    att   = 1-cos(el)^1.6
    att
end


