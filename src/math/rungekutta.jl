# Copyright (c) 2015-2017 Michael Eastwood
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

"The Butcher tableau for the 2nd-order Runge-Kutta method."
const RK2_tableau = @SMatrix [1/2 0.0
                              0.0 1.0]

"The Butcher tableau for the 3rd-order Runge-Kutta method."
const RK3_tableau = @SMatrix [1/2 0.0 0.0;
                             -1.0 2.0 0.0;
                              1/6 2/3 1/6]

"The Butcher tableau for the 4th-order Runge-Kutta method."
const RK4_tableau = @SMatrix [1/2 0.0 0.0 0.0;
                              0.0 1/2 0.0 0.0;
                              0.0 0.0 1.0 0.0;
                              1/6 1/3 1/3 1/6]

struct RKWorkspace{V<:AbstractVector}
    x :: V
    k :: Vector{V}
    δ :: V
end

function RKWorkspace(v::V, N::Int) where {V<:AbstractVector}
    x = similar(v)
    k = [similar(v) for n = 1:N]
    δ = similar(v)
    RKWorkspace(x, k, δ)
end

function rk2!(workspace, f!, x)
    k1 = workspace.k[1]; k2 = workspace.k[2]
    f!(k1, x)
    @. workspace.x = x + 0.5*k1
    f!(k2, workspace.x)
    @. workspace.δ = x + k2
    workspace.δ[:]
end

function rk4!(workspace, f!, x)
    k1 = workspace.k[1]; k2 = workspace.k[2]
    k3 = workspace.k[3]; k4 = workspace.k[4]
    f!(k1, x)
    @. workspace.x = x + 0.5*k1
    f!(k2, workspace.x)
    @. workspace.x = x + 0.5*k2
    f!(k3, workspace.x)
    @. workspace.x = x + k3
    f!(k4, workspace.x)
    @. workspace.δ = x + (k1 + 2*k2 + 2*k3 + k4)/6
    workspace.δ[:]
end

#"""
#    RK2(step)
#
#Wrap a step function with the 2nd-order Runge-Kutta method. For example
#
#    # Compute an approximation of exp(1)
#    rk = RK2(x -> x)
#    rk(1) + 1
#
#**Arguments:**
#
#* `step(x, args...)` must return the step to take from the starting location `x`
#* `x` is the starting location
#* `args` is simply passed on as the second argument to `step`
#"""
#@generate_step RK2 RK2_tableau
#
#"""
#    RK3(step)
#
#Wrap a step function with the 3rd-order Runge-Kutta method. For example
#
#    # Compute an approximation of exp(1)
#    rk = RK3(x -> x)
#    rk(1) + 1
#
#**Arguments:**
#
#* `step(x, args...)` must return the step to take from the starting location `x`
#* `x` is the starting location
#* `args` is simply passed on as the second argument to `step`
#"""
#@generate_step RK3 RK3_tableau
#
#"""
#    RK4(step)
#
#Wrap a step function with the 4th-order Runge-Kutta method. For example
#
#    # Compute an approximation of exp(1)
#    rk = RK4(x -> x)
#    rk(1) + 1
#
#**Arguments:**
#
#* `step(x, args...)` must return the step to take from the starting location `x`
#* `x` is the starting location
#* `args` is simply passed on as the second argument to `step`
#"""
#@generate_step RK4 RK4_tableau
#
#"""
#    iterate(step, maxiter, tolerance, x, args...)
#
#Repeatedly call `step(x, args...)` while updating the value of `x` until either the step size is
#sufficiently small or the maximum number of iterations is reached.
#"""
#function iterate(step, maxiter, tolerance, x, args...)
#    iter = 0
#    converged = false
#    while !converged && iter < maxiter
#        δ = step(x, args...)
#        if vecnorm(δ) < tolerance * vecnorm(x)
#            converged = true
#        end
#        for i in eachindex(x, δ)
#            x[i] += δ[i]
#        end
#        iter += 1
#    end
#    converged
#end

