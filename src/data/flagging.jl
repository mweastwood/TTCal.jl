# Copyright (c) 2015-2017 Michael Eastwood
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

function flag_antennas!(dataset::Dataset, antennas)
    for time = 1:Ntime(dataset), frequency = 1:Nfreq(dataset)
        visibilities = dataset[frequency, time]
        for antenna1 in antennas, antenna2 = 1:Nant(dataset)
            flag!(visibilities, antenna1, antenna2)
        end
    end
    antennas
end

function flag_baselines!(dataset::Dataset, baselines)
    for time = 1:Ntime(dataset), frequency = 1:Nfreq(dataset)
        visibilities = dataset[frequency, time]
        for (antenna1, antenna2) in baselines
            flag!(visibilities, antenna1, antenna2)
        end
    end
    baselines
end

function flag_frequencies!(dataset::Dataset, frequencies)
    for time = 1:Ntime(dataset), frequency in frequencies
        visibilities = dataset[frequency, time]
        for antenna1 = 1:Nant(dataset), antenna2 = antenna1:Nant(dataset)
            flag!(visibilities, antenna1, antenna2)
        end
    end
    frequencies
end

function match_flags!(to, from)
    for time = 1:Ntime(to), frequency = 1:Nfreq(to)
        to_visibilities   =   to[frequency, time]
        from_visibilities = from[frequency, time]
        for antenna1 = 1:Nant(to), antenna2 = antenna1:Nant(to)
            if isflagged(from_visibilities, antenna1, antenna2)
                flag!(to_visibilities, antenna1, antenna2)
            end
        end
    end
end

function flag_short_baselines!(dataset, minuvw)
    metadata = dataset.metadata
    for time = 1:Ntime(dataset), frequency = 1:Nfreq(dataset)
        visibilities = dataset[frequency, time]
        ν = metadata.frequencies[frequency]
        λ = u"c" / ν
        for antenna1 = 1:Nant(dataset), antenna2 = antenna1:Nant(dataset)
            baseline_vector = metadata.positions[antenna1] - metadata.positions[antenna2]
            baseline_length = norm(baseline_vector)
            if baseline_length < minuvw * λ
                flag!(visibilities, antenna1, antenna2)
            end
        end
    end
    dataset
end

