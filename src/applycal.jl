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

function applycal!(ms::Table,
                   gains::Array{Complex64};
                   apply_to_corrected::Bool = false,
                   force_imaging_columns::Bool = false)
    ant1,ant2 = ants(ms)
    data = apply_to_corrected? ms["CORRECTED_DATA"] : ms["DATA"]
    applycal!(data,gains,ant1,ant2)
    if force_imaging_columns || Tables.checkColumnExists(ms,"CORRECTED_DATA")
        ms["CORRECTED_DATA"] = data
    else
        ms["DATA"] = data
    end
    data
end

function applycal!(data::Array{Complex64,3},
                   gains::Array{Complex64,3},
                   ant1::Vector{Int32},
                   ant2::Vector{Int32})
    # gains are from bandpass(...)
    Nbase = length(ant1)
    Nfreq = size(gains,3)
    for α = 1:Nbase, β = 1:Nfreq
        data[1,β,α] = (gains[ant1[α],1,β]*conj(gains[ant2[α],1,β]))\data[1,β,α]
        data[2,β,α] = (gains[ant1[α],1,β]*conj(gains[ant2[α],2,β]))\data[2,β,α] # this one could be
        data[3,β,α] = (gains[ant1[α],2,β]*conj(gains[ant2[α],1,β]))\data[3,β,α] # swapped with this one
        data[4,β,α] = (gains[ant1[α],2,β]*conj(gains[ant2[α],2,β]))\data[4,β,α]
    end
end

function applycal!(data::Array{Complex64,3},
                   gains::Array{Complex64,4},
                   ant1::Vector{Int32},
                   ant2::Vector{Int32})
    # gains are from polcal(...)
    # TODO: speed this up (it is very slow)
    Nbase = length(ant1)
    Nfreq = size(gains,4)
    for α = 1:Nbase, β = 1:Nfreq
        V = [data[1,β,α] data[3,β,α]
             data[2,β,α] data[4,β,α]]
        G1 = gains[:,:,ant1[α],β]
        G2 = gains[:,:,ant2[α],β]
        V = G2'\(V/G1)
        data[1,β,α] = V[1,1]
        data[2,β,α] = V[2,1]
        data[3,β,α] = V[1,2]
        data[4,β,α] = V[2,2]
    end
end

