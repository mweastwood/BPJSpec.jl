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

struct Metadata
    lmax :: Int
    mmax :: Int
    ν :: Vector{typeof(1.0*u"Hz")}
    baselines :: Vector{Baseline}
    beam :: Map
    phase_center :: Direction
end

# compatibility with TTCal

function from_ttcal(frame, ttcal_metadata, lmax, mmax, ν, beam)
    baselines = ttcal_baselines(frame, ttcal_metadata)
    phase_center = ttcal_phase_center(frame, ttcal_metadata)
    Metadata(lmax, mmax, ν, baselines, beam, phase_center)
end

function ttcal_baselines(frame, ttcal_metadata)
    antennas = [measure(frame, antenna.position, pos"ITRF") for antenna in ttcal_metadata.antennas]
    baselines = Baseline[]
    for α = 1:length(ttcal_metadata.baselines)
        antenna1 = antennas[ttcal_metadata.baselines[α].antenna1]
        antenna2 = antennas[ttcal_metadata.baselines[α].antenna2]
        baseline = Baseline(baseline"ITRF",
                            antenna1.x - antenna2.x,
                            antenna1.y - antenna2.y,
                            antenna1.z - antenna2.z)
        push!(baselines, baseline)
    end
    baselines
end

function ttcal_phase_center(frame, ttcal_metadata)
    measure(frame, ttcal_metadata.phase_center, dir"ITRF")
end

