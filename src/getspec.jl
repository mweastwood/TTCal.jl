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
    getspec(ms::MeasurementSet, dir::Direction) -> xx,xy,yx,yy

This function extracts the spectrum in a given direction by means of an
inverse discrete Fourier transform.

Note that no gridding is performed,
so this does *not* use a fast Fourier transform. However, the inverse
discrete Fourier transform *is* the least squares estimator for the flux
in a given direction (if all baselines are weighted equally).
"""
function getspec(ms::MeasurementSet, dir::Direction)
    if !isabovehorizon(ms.frame,dir)
        error("Direction is below the horizon.")
    end

    l,m = dir2lm(ms.frame,ms.phase_direction,dir)
    data  = get_corrected_data(ms)
    flags = get_flags(ms)

    getspec(data,flags,l,m,
            ms.u,ms.v,ms.w,ms.ν,
            ms.ant1,ms.ant2)
end

doc"""
    getspec(data, flags, l, m, u, v, w, ν, ant1, ant2) -> xx,xy,yx,yy

Compute the spectrum of a source located at $(l,m)$ in all of the
polarized correlation products.
"""
function getspec(data::Array{Complex64,3},
                 flags::Array{Bool,3},
                 l,m,u,v,w,ν,ant1,ant2)
    fringe = fringepattern(l,m,u,v,w,ν)
    xx = zeros(length(ν))
    xy = zeros(length(ν))
    yx = zeros(length(ν))
    yy = zeros(length(ν))
    count  = zeros(Int,length(ν)) # The number of baselines used in the calculation
    for α = 1:length(u)
        # Don't use auto-correlations
        ant1[α] == ant2[α] && continue
        for β = 1:length(ν)
            any(slice(flags,:,β,α)) && continue
            # Taking the real part of A is equivalent to
            # computing (A + conj(A))/2. The conjugate of A,
            # in this case, is the baseline going in the
            # opposite direction.
            z = conj(fringe[β,α])
            xx[β] += real(data[1,β,α]*z)
            xy[β] += real(data[2,β,α]*z)
            yx[β] += real(data[3,β,α]*z)
            yy[β] += real(data[4,β,α]*z)
            count[β] += 1
        end
    end
    for β = 1:length(ν)
        xx[β] = xx[β]/count[β]
        xy[β] = xy[β]/count[β]
        yx[β] = yx[β]/count[β]
        yy[β] = yy[β]/count[β]
    end
    xx,xy,yx,yy
end

