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

function fringepattern(l,m,u,v,w,ν)
    output = Array(Complex64,length(ν),length(u))
    fringepattern!(output,l,m,u,v,w,ν)
    output
end

function fringepattern!{T<:Complex}(output::Array{T,2},l,m,u,v,w,ν)
    n = sqrt(1-l^2-m^2)-1
    Δν = ν[2] - ν[1]
    Nbase = length(u)
    Nfreq = length(ν)

    fringe = Array(Complex64,Nfreq)
    for α = 1:Nbase
        # Get the fringe pattern for the baseline
        τ = 2π*(u[α]*l+v[α]*m+w[α]*n)/c
        ϕ = τ*ν[1]
        Δϕ = τ*Δν
        fringepattern!(fringe,ϕ,Δϕ)
        for β = 1:Nfreq
            output[β,α] = fringe[β]
        end
    end
    nothing
end

doc"""
    fringepattern(ϕ, Δϕ, N)

Compute $\exp(i(\phi+n\Delta\phi))$ where $\phi$, $\Delta\phi$, and $n = 1,\ldots,N$ define
an equally spaced grid of points.

Using the sine and cosine angle addition rules, you can define
an iterative method such that you only need to compute sines
and cosines for a single iteration.
"""
function fringepattern(ϕ,Δϕ,N::Int)
    output = Array(Complex64,N)
    fringepattern!(output,ϕ,Δϕ)
    output
end

function fringepattern!{T<:Complex}(output::Array{T,1},ϕ,Δϕ)
    N = length(output)
    sin_Δϕ = sin(Δϕ)
    cos_Δϕ = cos(Δϕ)
    output[1] = complex(cos(ϕ),sin(ϕ))
    for n = 1:N-1
        output[n+1] = complex(real(output[n])*cos_Δϕ - imag(output[n])*sin_Δϕ,
                              imag(output[n])*cos_Δϕ + real(output[n])*sin_Δϕ)
    end
    output
end

