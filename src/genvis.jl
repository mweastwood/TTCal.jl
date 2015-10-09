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
    genvis(ms::Table,
           sources::Vector{PointSource},
           beam::BeamModel)

Generate model visibilities for the given list of sources
and the given beam model.

No gridding is performed, so the runtime of this naive
algorithm scales as $O(N_{base} \times N_{source})$.
"""
function genvis(ms::Table,
                sources::Vector{PointSource},
                beam::BeamModel)
    phase_dir = MeasurementSets.phase_direction(ms)
    u,v,w = MeasurementSets.uvw(ms)
    ν = MeasurementSets.frequency(ms)

    frame = ReferenceFrame()
    set!(frame,MeasurementSets.position(ms))
    set!(frame,MeasurementSets.time(ms))

    genvis(frame,phase_dir,sources,beam,u,v,w,ν)
end

function genvis(frame::ReferenceFrame,
                phase_dir::Direction,
                sources::Vector{PointSource},
                beam::BeamModel,
                u,v,w,ν)
    model = zeros(Complex64,4,length(ν),length(u))
    for source in sources
        genvis!(model,frame,phase_dir,source,beam,u,v,w,ν)
    end
    model
end

function genvis!(model::Array{Complex64,3},
                 frame::ReferenceFrame,
                 phase_dir::Direction,
                 source::PointSource,
                 beam::BeamModel,
                 u,v,w,ν)
    az,el  = azel(frame,source)
    l,m    = lm(frame,phase_dir,source)
    fringe = fringepattern(l,m,u,v,w,ν)

    # Get the Stokes parameters and apply the beam model
    correlations = zeros(Complex128,4,length(ν))
    for β = 1:length(ν)
        M = mueller(beam(ν,az,el))
        IQUV = stokes_flux(source,ν[β])
        correlations[:,β] = linear(M*IQUV)
    end

    genvis!(model,fringe,correlations)
end

function genvis!(model::Array{Complex64,3},
                 fringe,correlations)
    Nfreq,Nbase = size(fringe)
    @inbounds for α = 1:Nbase, β = 1:Nfreq
        model[1,β,α] += correlations[1,β]*fringe[β,α] # xx
        model[2,β,α] += correlations[2,β]*fringe[β,α] # xy
        model[3,β,α] += correlations[3,β]*fringe[β,α] # yx
        model[4,β,α] += correlations[4,β]*fringe[β,α] # yy
    end
    model
end

