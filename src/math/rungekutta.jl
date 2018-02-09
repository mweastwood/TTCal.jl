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

