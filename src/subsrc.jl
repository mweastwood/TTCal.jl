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

function subsrc!(ms::Table,sources::Vector{PointSource})
    frame = reference_frame(ms)
    u,v,w = uvw(ms)
    ν     = freq(ms)
    data  = Tables.checkColumnExists(ms,"CORRECTED_DATA")? ms["CORRECTED_DATA"] : ms["DATA"]

    subtracted = subsrc(frame,data,u,v,w,ν,sources)
    ms["CORRECTED_DATA"] = subtracted
    subtracted
end

################################################################################
# Internal Interface

function subsrc(frame::ReferenceFrame,
                data::Array{Complex64,3},
                u::Vector{quantity(Float64,Meter)},
                v::Vector{quantity(Float64,Meter)},
                w::Vector{quantity(Float64,Meter)},
                ν::Vector{quantity(Float64,Hertz)},
                sources::Vector{PointSource})
    sources = filter(sources) do source
        isabovehorizon(frame,source)
    end

    model = genvis(frame,sources,u,v,w,ν)
    # Re-use the space allocated for the model visibilities
    # to store the model subtracted visibilities.
    for i = 1:length(model)
        model[i] = data[i]-model[i]
    end
    model
end

