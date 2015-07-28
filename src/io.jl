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

function write_gains(filename,gains::Array{Complex64,3},gain_flags::Array{Bool,3})
    # gains are from bandpass(...)
    Nant,Npol,Nchan = size(gains)
    open(filename,"w") do f
        write(f,'B')
        write(f,Int32(Nant))
        write(f,Int32(Nchan))
        write(f,permutedims(gain_flags,(2,3,1)))
        write(f,permutedims(gains,(2,3,1)))
    end
end

function write_gains(filename,gains::Array{Complex128,4},gain_flags::Array{Bool,2})
    # gains are from polcal(...)
    Npol1,Npol2,Nant,Nchan = size(gains)
    open(filename,"w") do f
        write(f,'J')
        write(f,Int32(Nant))
        write(f,Int32(Nchan))
        write(f,permutedims(gain_flags,(2,1)))
        write(f,permutedims(gains,(1,2,4,3)))
    end
end

function read_gains(filename)
    local gains,gain_flags
    open(filename,"r") do f
        T = read(f,Char)
        Nant = Int(read(f,Int32))
        Nchan = Int(read(f,Int32))
        if T == 'B'
            gain_flags = permutedims(read(f,Bool,(2,Nchan,Nant)),(3,1,2))
            gains = permutedims(read(f,Complex64,(2,Nchan,Nant)),(3,1,2))
        elseif T == 'J'
            gain_flags = permutedims(read(f,Bool,(Nchan,Nant)),(2,1))
            gains = permutedims(read(f,Complex128,(2,2,Nchan,Nant)),(1,2,4,3))
        else
            error("Unknown calibration type.")
        end
    end
    gains,gain_flags
end

