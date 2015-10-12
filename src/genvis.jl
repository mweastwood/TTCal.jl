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

doc"""
    genvis(ms::MeasurementSet,
           sources::Union{PointSource,Vector{PointSource}},
           beam::BeamModel)

Generate model visibilities for the given list of sources
and the given beam model.

No gridding is performed, so the runtime of this naive
algorithm scales as $O(N_{base} \times N_{source})$.
"""
function genvis(ms::MeasurementSet,
                sources::Union{PointSource,Vector{PointSource}},
                beam::BeamModel)
    genvis(ms.frame,ms.phase_direction,
           sources,beam,ms.u,ms.v,ms.w,ms.ν)
end

function genvis(frame::ReferenceFrame,
                phase_direction::Direction,
                sources::Vector{PointSource},
                beam::BeamModel,
                u,v,w,ν)
    model = zeros(Complex64,4,length(ν),length(u))
    for source in sources
        genvis!(model,frame,phase_direction,source,beam,u,v,w,ν)
    end
    model
end

function genvis(frame::ReferenceFrame,
                phase_direction::Direction,
                source::PointSource,
                beam::BeamModel,
                u,v,w,ν)
    model = zeros(Complex64,4,length(ν),length(u))
    genvis!(model,frame,phase_direction,source,beam,u,v,w,ν)
    model
end

function genvis!(model::Array{Complex64,3},
                 frame::ReferenceFrame,
                 phase_direction::Direction,
                 source::PointSource,
                 beam::BeamModel,
                 u,v,w,ν)
    az,el  = azel(frame,source)
    l,m    = lm(frame,phase_direction,source)
    fringe = fringepattern(l,m,u,v,w,ν)

    # Get the Stokes parameters and apply the beam model
    xx = zeros(length(ν))
    xy = zeros(length(ν))
    yx = zeros(length(ν))
    yy = zeros(length(ν))
    for β = 1:length(ν)
        M = mueller(beam(ν[β],az,el))
        IQUV = stokes_flux(source,ν[β])
        correlations = linear(M*IQUV)
        xx[β] = correlations[1]
        xy[β] = correlations[2]
        yx[β] = correlations[3]
        yy[β] = correlations[4]
    end

    genvis!(model,xx,xy,yx,yy,l,m,u,v,w,ν)
end

function genvis(xx,xy,yx,yy,
                l,m,u,v,w,ν)
    model = zeros(Complex64,4,length(ν),length(u))
    genvis!(model,xx,xy,yx,yy,l,m,u,v,w,ν)
    model
end

function genvis!(model::Array{Complex64,3},
                 xx,xy,yx,yy,
                 l,m,u,v,w,ν)
    fringe = fringepattern(l,m,u,v,w,ν)
    Nfreq,Nbase = size(fringe)
    @inbounds for α = 1:Nbase, β = 1:Nfreq
        model[1,β,α] += xx[β]*fringe[β,α]
        model[2,β,α] += xy[β]*fringe[β,α]
        model[3,β,α] += yx[β]*fringe[β,α]
        model[4,β,α] += yy[β]*fringe[β,α]
    end
    model
end

