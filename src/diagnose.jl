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

@doc """
Diagnose a calibration by looking for bad antennas.

The algorithms employed here to identify misbehaving antennas are
tailored for the OVRO LWA and may not be ideal for other interferometers.
""" ->
function diagnose(filename)
    gains,gain_flags = read_gains(filename)
    diagnose(gains,gain_flags)
end

function diagnose(gains::Array{Complex64,3},
                  gain_flags::Array{Bool,3};
                  Nprint=20)
    Nant = size(gains,1)

    # Calculate the correlation function of the antenna gains
    xcorr = zeros(Nant,Nant)
    ycorr = zeros(Nant,Nant)
    for ant1 = 1:Nant, ant2 = 1:Nant
        x1 = abs(slice(gains,ant1,1,:))
        x2 = abs(slice(gains,ant2,1,:))
        y1 = abs(slice(gains,ant1,2,:))
        y2 = abs(slice(gains,ant2,2,:))
        σx1 = mean((x1-mean(x1)) .* conj(x1-mean(x1)))
        σx2 = mean((x2-mean(x2)) .* conj(x2-mean(x2)))
        σy1 = mean((y1-mean(y1)) .* conj(y1-mean(y1)))
        σy2 = mean((y2-mean(y2)) .* conj(y2-mean(y2)))
        σx1x2 = mean((x1-mean(x1)) .* conj(x2-mean(x2)))
        σy1y2 = mean((y1-mean(y1)) .* conj(y2-mean(y2)))
        xcorr[ant1,ant2] = σx1x2/sqrt(σx1*σx2)
        ycorr[ant1,ant2] = σy1y2/sqrt(σy1*σy2)
    end
    xcorr[isnan(xcorr)] = 0.0
    ycorr[isnan(ycorr)] = 0.0

    # Factor the correlation matriices
    λ,v = eigs(xcorr,nev=1,which=:LM)
    x = sqrt(λ[1])*v[:,1]
    λ,v = eigs(ycorr,nev=1,which=:LM)
    y = sqrt(λ[1])*v[:,1]

    # Try to make all the good antennas positive
    # (this assumes that >50% of the antennas are good)
    x *= sign(median(x))
    y *= sign(median(y))
    scores = [x;y]

    # Print the worst antennas
    println("Consider flagging one or more of the following antennas:")
    println(" Antenna | Pol | Score ")
    println("---------|-----|-------")
    badants = sortperm(scores)
    N = 1
    idx = 1
    while N < Nprint
        ant = badants[idx]
        pol = 1
        if ant > Nant
            ant -= Nant
            pol = 2
        end
        if all(gain_flags[ant,pol,:])
            idx += 1
            continue
        end
        @printf("   %3d   |  %s  | %+4.2f \n",ant,pol == 1?"a":"b",scores[badants[idx]])
        N += 1
        idx += 1
    end
end

