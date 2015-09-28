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

####################################################################
# Most basic of beam models
"""
Generate beam model to attenuate source model flux (to be used in genvis.jl).
As an initial attempt, uses the naive sin(el)^1.6 beam, normalized to 1 at zenith.
"""
function beammodel(frame,phase_dir,l,m)
    dir   = lm2dir(phase_dir,l,m)
    az,el = dir2azel(frame,dir)

    att   = sin(el)^1.6
    att
end


####################################################################
# J. Dowell empirical beam model, function of frequency
"""
Generate beam model to attenuate source model flux (to be used in genvis.jl).
Using the J. Dowell beam model (http://www.faculty.ece.vt.edu/swe/lwa/memo/lwa0178a.pdf).
"""
file    = jldopen("/scr/mmanders/LWA/ttcal_tutorial/lwa1-dipole-emp.jld", "r")
jdmodel = JLD.read(file, "jdmodel")
fitX    = jdmodel["fitX"]
fitY    = jdmodel["fitY"]

function beammodel(frame,phase_dir,l,m,ν)
    dir   = lm2dir(phase_dir,l,m)
    az,el = dir2azel(frame,dir)
    # zenith angle
    ze    = π/2-el

    px = sqrt( ( planeE(fitX,ze,ν) * cos(az) )^2 + ( planeH(fitX,ze,ν) * sin(az) )^2 )
    py = sqrt( ( planeE(fitY,ze,ν) * cos(az) )^2 + ( planeH(fitY,ze,ν) * sin(az) )^2 )
    p  = (px + py) / 2
    p
end

function planeE(fitPol,ze,ν)
    alpha = polyval(fitPol[1,1,:],ν)
    beta  = polyval(fitPol[1,2,:],ν)
    gamma = polyval(fitPol[1,3,:],ν)
    delta = polyval(fitPol[1,4,:],ν)

    out = (1 - (2*ze/π)^alpha) * cos(ze)^beta
    out += gamma * (2*ze/π) * cos(ze)^delta
    out
end

function planeH(fitPol,ze,ν)
    alpha = polyval(fitPol[2,1,:],ν)
    beta  = polyval(fitPol[2,2,:],ν)
    gamma = polyval(fitPol[2,3,:],ν)
    delta = polyval(fitPol[2,4,:],ν)
    
    out = (1 - (2*ze/π)^alpha) * cos(ze)^beta
    out += gamma * (2*ze/π) * cos(ze)^delta
    out
end


"""
Evaluate polynomial of degree N-1 where N = len(coeffs).
"""
function polyval(coeffs,param)
    y = 0
    for c = coeffs
        y = y * param + c
    end
    y
end
