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
    genvis(ms::MeasurementSet, sources, beam::BeamModel)

Generate model visibilities for the given list of sources
and the given beam model.

No gridding is performed, so the runtime of this naive
algorithm scales as $O(N_{base} \times N_{source})$.
"""
function genvis(ms::MeasurementSet, sources, beam::BeamModel)
    model = zeros(Complex64,4,ms.Nfreq,ms.Nbase)
    genvis!(model,ms.frame,ms.phase_direction,sources,beam,ms.u,ms.v,ms.w,ms.ν)
    model
end

function genvis!(model::Array{Complex64,3},
                 frame::ReferenceFrame,
                 phase_direction::Direction,
                 sources::Vector{Source},
                 beam::BeamModel,
                 u,v,w,ν)
    for source in sources
        genvis!(model,frame,phase_direction,source,beam,u,v,w,ν)
    end
    model
end

function genvis!(model::Array{Complex64,3},
                 frame::ReferenceFrame,
                 phase_direction::Direction,
                 source::Source,
                 beam::BeamModel,
                 u,v,w,ν)
    for component in source.components
        genvis!(model,frame,phase_direction,component,beam,u,v,w,ν)
    end
    model
end

function genvis!(model::Array{Complex64,3},
                 frame::ReferenceFrame,
                 phase_direction::Direction,
                 point::Point,
                 beam::BeamModel,
                 u,v,w,ν)
    Nfreq  = length(ν)
    az,el  = local_azel(frame,point)
    l,m    = direction_cosines(phase_direction,measure(frame,point.direction,dir"J2000"))
    fringe = fringepattern(l,m,u,v,w,ν)

    flux = Array{HermitianJonesMatrix}(Nfreq)
    for β = 1:length(ν)
        jones  = beam(ν[β],az,el)
        stokes = point.spectrum(ν[β])
        uncorrupted = linear(stokes) # convert from I,Q,U,V to xx,xy,yx,yy
        corrupted   = congruence_transform(jones,uncorrupted)
        flux[β] = corrupted
    end

    genvis!(model,flux,l,m,u,v,w,ν)
end

function genvis(flux,l,m,u,v,w,ν)
    model = zeros(Complex64,4,length(ν),length(u))
    genvis!(model,flux,l,m,u,v,w,ν)
    model
end

function genvis!(model::Array{Complex64,3},
                 flux,l,m,u,v,w,ν)
    fringe = fringepattern(l,m,u,v,w,ν)
    Nfreq,Nbase = size(fringe)
    @inbounds for α = 1:Nbase, β = 1:Nfreq
        model[1,β,α] += flux[β].xx*fringe[β,α]
        model[2,β,α] += flux[β].xy*fringe[β,α]
        model[3,β,α] += conj(flux[β].xy)*fringe[β,α]
        model[4,β,α] += flux[β].yy*fringe[β,α]
    end
    model
end

