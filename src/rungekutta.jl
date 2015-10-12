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

abstract StepFunction

"""
The Butcher tableau for the 2nd-order Runge-Kutta method.
"""
const RK2_tableau = [1/2 0.0
                     0.0 1.0]

"""
The Butcher tableau for the 3rd-order Runge-Kutta method.
"""
const RK3_tableau = [1/2 0.0 0.0;
                    -1.0 2.0 0.0;
                     1/6 2/3 1/6]

"""
The Butcher tableau for the 4th-order Runge-Kutta method.
"""
const RK4_tableau = [1/2 0.0 0.0 0.0;
                     0.0 1/2 0.0 0.0;
                     0.0 0.0 1.0 0.0;
                     1/6 1/3 1/3 1/6]

"""
This singleton type is used to indicate which Runge-Kutta method
should be used. For example, `RK{4}` tells us to use the RK4 method.
"""
immutable RK{N} end
RK(N) = RK{N}()
rkorder{N}(::RK{N}) = N
rkorder{N}(::Type{RK{N}}) = N

const RK2 = RK(2)
const RK3 = RK(3)
const RK4 = RK(4)

@generated function rkstep(step::StepFunction,rk::RK,x,args...)
    N = rkorder(rk)
    T = eltype(x)
    Ndims = ndims(x)
    tableau = symbol("RK$(N)_tableau")

    quote
        # Calculate the intermediate steps k
        k = Array{Array{$T,$Ndims}}($N)
        k[1] = step(x,args...)
        for row = 1:$N-1
            x′ = copy(x)
            for col = 1:row
                $tableau[row,col] == 0 && continue
                k′ = k[col]
                for i in eachindex(x′)
                    @inbounds x′[i] += k′[i]*$tableau[row,col]
                end
            end
            k[row+1] = step(x′,args...)
        end
        # Calculate the final step δ
        δ = zeros(eltype(x),size(x))
        for col = 1:$N
            k′ = k[col]
            for i in eachindex(δ)
                @inbounds δ[i] += k′[i]*$tableau[$N,col]
            end
        end
        δ
    end
end

@doc """
Take a Runge-Kutta step.
* `step(x,args...)` must return the list of steps to take from the given location `x`
* `rk` is the order of the Runge-Kutta method to use (eg. `RK4`)
* `x` is the starting location
* `args` is simply passed on as the second argument to `step`
""" rkstep

macro iterate(step,rk,maxiter,tolerance,x,args...)
    quote
        iter = 0
        converged = false
        while !converged && iter < $maxiter
            δ = rkstep($step,$rk,$x,$(args...))
            if vecnorm(δ) < $tolerance*vecnorm($x)
                converged = true
            end
            for i in eachindex($x)
                $x[i] += δ[i]
            end
            iter += 1
        end
        converged
    end
end

