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

################################################################################
# Public Interface

@doc """
This function extracts the spectrum in a given direction by means of an
inverse discrete Fourier transform. Note that no gridding is performed,
so this does *not* use a fast Fourier transform. However, the inverse
discrete Fourier transform *is* the least squares estimator for the flux
in a given direction (if all baselines are weighted equally).
""" ->
function getspec(ms::Table,
                 dir::Direction)
    data  = ms["CORRECTED_DATA"]
    flags = ms["FLAG"]
    frame = reference_frame(ms)
    l,m = dir2lm(frame,dir)
    u,v,w = uvw(ms)
    ν = freq(ms)
    ant1,ant2 = ants(ms)
    getspec_internal(data,flags,l,m,u,v,w,ν,ant1,ant2)
end

getspec(ms::Table,source::PointSource) = getspec(ms,direction(source))

################################################################################
# Internal Interface

function getspec_internal(data::Array{Complex64,3},
                          flags::Array{Bool,3},
                          l::Float64,
                          m::Float64,
                          u::Vector{quantity(Float64,Meter)},
                          v::Vector{quantity(Float64,Meter)},
                          w::Vector{quantity(Float64,Meter)},
                          ν::Vector{quantity(Float64,Hertz)},
                          ant1::Vector{Int32},
                          ant2::Vector{Int32})
    fringe = fringepattern(l,m,u,v,w,ν)
    xx = zeros(Float64,length(ν))
    xy = zeros(Complex128,length(ν))
    yy = zeros(Float64,length(ν))
    count  = zeros(Int,length(ν)) # The number of baselines used in the calculation
    for α = 1:length(u)
        # Don't use auto-correlations
        ant1[α] == ant2[α] && continue
        for β = 1:length(ν)
            any(flags[:,β,α]) && continue
            # Taking the real part of A is equivalent to
            # computing (A + conj(A))/2. The conjugate of A,
            # in this case, is the baseline going in the
            # opposite direction. Including this information
            # constrains the spectrum to be real.
            z = conj(fringe[β,α])
            xx[β] += real(data[1,β,α]*z) # xx
            xy[β] += 0.5*(data[2,β,α]*z + conj(data[3,β,α]*z)) # xy
            yy[β] += real(data[4,β,α]*z) # yy
            count[β] += 1
        end
    end
    xx = xx ./ count
    xy = xy ./ count
    yy = yy ./ count
    xx,xy,yy
end

function xy2stokes(xx,xy,yy)
    I = 0.5*(xx+yy)
    Q = 0.5*(yy-xx)
    U = real(xy)
    V = -imag(xy)
    I,Q,U,V
end

