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
    getspec(ms::MeasurementSet, direction::Direction;
            minuvw = 0.0) -> Vector{HermitianJonesMatrix}

This function extracts the spectrum in a given direction by means of an
inverse discrete Fourier transform.

Note that no gridding is performed,
so this does *not* use a fast Fourier transform. However, the inverse
discrete Fourier transform *is* the least squares estimator for the flux
in a given direction (if all baselines are weighted equally).
"""
function getspec(ms::MeasurementSet, direction::Direction;
                 minuvw::Float64 = 0.0)
    if !isabovehorizon(ms.frame,direction)
        error("Direction is below the horizon.")
    end

    j2000 = measure(ms.frame,direction,dir"J2000")
    l,m   = direction_cosines(ms.phase_direction,j2000)

    data  = get_corrected_data(ms)
    flags = get_flags(ms)
    flag_short_baselines!(flags,minuvw,ms.u,ms.v,ms.w,ms.ν)

    getspec(data,flags,l,m,
            ms.u,ms.v,ms.w,ms.ν,
            ms.ant1,ms.ant2)
end

doc"""
    getspec(data, flags, l, m, u, v, w, ν, ant1, ant2) -> Vector{HermitianJonesMatrix}

Compute the spectrum of a source located at $(l,m)$ in all of the
polarized correlation products.
"""
function getspec(data::Array{Complex64,3},
                 flags::Array{Bool,3},
                 l,m,u,v,w,ν,ant1,ant2)
    Nfreq  = length(ν)
    fringe = fringepattern(l,m,u,v,w,ν)
    flux   = zeros(HermitianJonesMatrix,Nfreq)
    count  = zeros(Int,Nfreq) # The number of baselines used in the calculation
    for α = 1:length(u)
        # Don't use auto-correlations
        ant1[α] == ant2[α] && continue
        for β = 1:length(ν)
            any(slice(flags,:,β,α)) && continue
            # Taking the real part of A is equivalent to
            # computing (A + conj(A))/2. The conjugate of A,
            # in this case, is the baseline going in the
            # opposite direction.
            f = conj(fringe[β,α])
            F = HermitianJonesMatrix(real(data[1,β,α]*f), # xx
                                     complex(real(data[2,β,α]*f),-imag(data[3,β,α]*f)), #xy
                                     real(data[4,β,α]*f)) # yy
            flux[β]  += F
            count[β] += 1
        end
    end
    for β = 1:length(ν)
        flux[β] /= count[β]
    end
    flux
end

