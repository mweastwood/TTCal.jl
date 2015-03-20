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
function genvis(ms::Table,sources::Vector{PointSource})
    u,v,w = uvw(ms)
    ν = freq(ms)
    frame = reference_frame(ms)
    genvis(frame,sources,u,v,w,ν)
end

################################################################################
# Internal Interface

function genvis(frame::ReferenceFrame,
                sources::Vector{PointSource},
                u::Vector{quantity(Float64,Meter)},
                v::Vector{quantity(Float64,Meter)},
                w::Vector{quantity(Float64,Meter)},
                ν::Vector{quantity(Float64,Hertz)})
    model = zeros(Complex64,4,length(ν),length(u))
    genvis!(model,frame,sources,u,v,w,ν)
    model
end

function genvis!(model::Array{Complex64,3},
                 frame::ReferenceFrame,
                 sources::Vector{PointSource},
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
                 source::PointSource,
                 u::Vector{quantity(Float64,Meter)},
                 v::Vector{quantity(Float64,Meter)},
                 w::Vector{quantity(Float64,Meter)},
                 ν::Vector{quantity(Float64,Hertz)})
    l,m = getlm(frame,source)
    I = getstokesI(source,ν)
    Q = getstokesQ(source,ν)
    U = getstokesU(source,ν)
    V = getstokesV(source,ν)
    genvis!(model,I,Q,U,V,l,m,u,v,w,ν)
end

function genvis!(model::Array{Complex64,3},
                 I::Vector{Float64},
                 Q::Vector{Float64},
                 U::Vector{Float64},
                 V::Vector{Float64},
                 l::Float64,
                 m::Float64,
                 u::Vector{quantity(Float64,Meter)},
                 v::Vector{quantity(Float64,Meter)},
                 w::Vector{quantity(Float64,Meter)},
                 ν::Vector{quantity(Float64,Hertz)})
    fringe = fringepattern(l,m,u,v,w,ν)
    for α = 1:length(u), β = 1:length(ν)
        # Based on experiments imaging model visibilities with wsclean
        # it appears that in order to have the model flux be equal to
        # the flux in the image, there is no factor of 0.5 in the
        # following statements.
        #
        # The implication of this is that you need the factor of 0.5
        # in the definition of the Stokes parameters such that, for
        # example I = 0.5*(xx+yy)
        model[1,β,α] += (I[β]-Q[β])*fringe[β,α] # xx
        model[2,β,α] += (U[β]-1im*V[β])*fringe[β,α] # xy
        model[3,β,α] += (U[β]+1im*V[β])*fringe[β,α] # yx
        model[4,β,α] += (I[β]+Q[β])*fringe[β,α] # yy
    end
    model
end

