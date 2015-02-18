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
Generate model visibilities for a given point source model.
No gridding is performed, so the runtime of this naive
algorithm scales as O(Nbase×Nsource).
""" ->
function genvis{T<:AbstractSource}(ms::Table,sources::Vector{T})
    u,v,w = uvw(ms)
    ν = freq(ms)
    frame = reference_frame(ms)
    genvis(frame,sources,u,v,w,ν)
end

################################################################################
# Internal Interface

function genvis{T<:AbstractSource}(frame::ReferenceFrame,
                                   sources::Vector{T},
                                   u::Vector{quantity(Float64,Meter)},
                                   v::Vector{quantity(Float64,Meter)},
                                   w::Vector{quantity(Float64,Meter)},
                                   ν::Vector{quantity(Float64,Hertz)})
    model = zeros(Complex64,4,length(ν),length(u))
    genvis!(model,frame,sources,u,v,w,ν)
    model
end

function genvis!{T<:AbstractSource}(model::Array{Complex64,3},
                                    frame::ReferenceFrame,
                                    sources::Vector{T},
                                    u::Vector{quantity(Float64,Meter)},
                                    v::Vector{quantity(Float64,Meter)},
                                    w::Vector{quantity(Float64,Meter)},
                                    ν::Vector{quantity(Float64,Hertz)})
    for source in sources
        genvis!(model,frame,source,u,v,w,ν)
    end
    model
end

function genvis!(model::Array{Complex64,3},
                 frame::ReferenceFrame,
                 source::AbstractSource,
                 u::Vector{quantity(Float64,Meter)},
                 v::Vector{quantity(Float64,Meter)},
                 w::Vector{quantity(Float64,Meter)},
                 ν::Vector{quantity(Float64,Hertz)})
    l,m = lm(frame,source)
    S = flux(source,ν) * sqrt(1-l^2-m^2)^1.9
    genvis!(model,S,l,m,u,v,w,ν)
end

function genvis!(model::Array{Complex64,3},
                 flux::Vector{Float64},
                 l::Float64,
                 m::Float64,
                 u::Vector{quantity(Float64,Meter)},
                 v::Vector{quantity(Float64,Meter)},
                 w::Vector{quantity(Float64,Meter)},
                 ν::Vector{quantity(Float64,Hertz)})
    fringe = fringepattern(l,m,u,v,w,ν)
    for α = 1:length(u), β = 1:length(ν)
        model[1,β,α] += 0.5*flux[β]*fringe[β,α] # xx
        model[4,β,α] += 0.5*flux[β]*fringe[β,α] # yy
    end
    model
end

