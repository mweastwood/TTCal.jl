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

using PyPlot

function peel(ms::Table,
              sources::Vector{PointSource},
              criteria::StoppingCriteria)
    frame = reference_frame(ms)
    dir   = phase_dir(ms)
    u,v,w = uvw(ms)
    ν = freq(ms)
    ant1,ant2 = ants(ms)
    data = corrected_data(ms)
    data_flags = flags(ms)

    Nsource = length(sources)
    Nfreq   = length(ν)
    Nant    = numrows(Table(ms[kw"ANTENNA"]))
    Nbase   = length(u)

    sources = filter(source -> isabovehorizon(frame,source),sources)
    coherencies = [genvis(dir,source,u,v,w,ν) for source in sources]
    #gains = [identity_matrices(Nant,Nfreq) for source in sources]
    #gain_flags = zeros(Bool,Nant,Nfreq)
    gains = [ones(Complex64,Nant,2,Nfreq) for source in sources]
    gain_flags = zeros(Bool,Nant,2,Nfreq)

    # Subtract all of the sources
    # (assuming the beam is unity towards each source)
    for coherency in coherencies
        for i in eachindex(data)
            data[i] -= coherency[i]
        end
    end

    # Derive a calibration towards each source
    for iter = 1:3
        @time for s = 1:Nsource
            tic()
            @show sources[s]
            coherency = coherencies[s]
            gains_toward_source = gains[s]

            corrupted = copy(coherency)
            corrupt!(corrupted,gains_toward_source,ant1,ant2)
            for i in eachindex(data)
                data[i] += corrupted[i]
            end

            #polcal!(gains_toward_source,gain_flags,
                    #data,coherency,data_flags,
                    #ant1,ant2,criteria,1/1000)
            bandpass!(gains_toward_source,gain_flags,
                      data,coherency,data_flags,
                      ant1,ant2,criteria)
            @show gains_toward_source[1]

            # Take the source back out of the measured visibilities,
            # but this time subtract it with the gains toward that direction.
            corrupted = copy(coherency)
            corrupt!(corrupted,gains_toward_source,ant1,ant2)
            for i in eachindex(data)
                data[i] -= corrupted[i]
            end
            toc()
        end
    end
    #applycal!(data,data_flags,gains[2],gain_flags,ant1,ant2)

    ms["CORRECTED_DATA"] = data
end

