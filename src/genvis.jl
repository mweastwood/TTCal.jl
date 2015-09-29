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

"""
Generate model visibilities for a given point source model.
No gridding is performed, so the runtime of this naive
algorithm scales as O(Nbase×Nsource).
"""
function genvis(ms::Table,sources::Vector{PointSource})
    phase_dir = MeasurementSets.phase_direction(ms)
    u,v,w = MeasurementSets.uvw(ms)
    ν = MeasurementSets.frequency(ms)

    frame = ReferenceFrame()
    set!(frame,MeasurementSets.position(ms))
    set!(frame,MeasurementSets.time(ms))

    genvis(frame,phase_dir,sources,u,v,w,ν)
end

genvis(ms::Table,source::PointSource) = genvis(ms,[source])

################################################################################
# Internal Interface

genvis(phase_dir::Direction,source::PointSource,u,v,w,ν) = genvis(phase_dir,[source],u,v,w,ν)

function genvis(frame::ReferenceFrame,
                phase_dir::Direction,
                sources::Vector{PointSource},
                u::Vector{Float64},
                v::Vector{Float64},
                w::Vector{Float64},
                ν::Vector{Float64})
    model = zeros(Complex64,4,length(ν),length(u))
    genvis!(model,frame,phase_dir,sources,u,v,w,ν)
    model
end

function genvis!(model::Array{Complex64,3},
                 frame::ReferenceFrame,
                 phase_dir::Direction,
                 sources::Vector{PointSource},
                 u::Vector{Float64},
                 v::Vector{Float64},
                 w::Vector{Float64},
                 ν::Vector{Float64})
    for source in sources
        genvis!(model,frame,phase_dir,source,u,v,w,ν)
    end
    model
end

function genvis!(model::Array{Complex64,3},
                 frame::ReferenceFrame,
                 phase_dir::Direction,
                 source::PointSource,
                 u::Vector{Float64},
                 v::Vector{Float64},
                 w::Vector{Float64},
                 ν::Vector{Float64})
    l,m = lm(frame,phase_dir,source)
    # att = beammodel(frame,phase_dir,l,m,mean(ν)) # currently using J.Dowell beam model
    att = beammodel(frame,phase_dir,l,m)
    I = stokesI(source,ν) * att
    Q = stokesQ(source,ν) * att
    U = stokesU(source,ν) * att
    V = stokesV(source,ν) * att
    genvis!(model,I,Q,U,V,l,m,u,v,w,ν)
end

function genvis!(model::Array{Complex64,3},
                 I::Vector{Float64},
                 Q::Vector{Float64},
                 U::Vector{Float64},
                 V::Vector{Float64},
                 l::Float64,
                 m::Float64,
                 u::Vector{Float64},
                 v::Vector{Float64},
                 w::Vector{Float64},
                 ν::Vector{Float64})
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

