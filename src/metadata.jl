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
    frequencies  :: Vector{typeof(1.0*u"Hz")}
    position     :: Position         # ITRF
    baselines    :: Vector{Baseline} # ITRF
    phase_center :: Direction        # ITRF
end

# compatibility with TTCal

function from_ttcal(ttcal_metadata)
    frequencies  = ttcal_metadata.frequencies
    position     = ttcal_position(ttcal_metadata)
    baselines    = ttcal_baselines(ttcal_metadata)
    phase_center = ttcal_metadata.phase_centers[1]
    Metadata(frequencies, position, baselines, phase_center)
end

function ttcal_position(ttcal_metadata)
    mean(ttcal_metadata.positions)
end

function ttcal_baselines(ttcal_metadata)
    positions = ttcal_metadata.positions
    baselines = Baseline[]
    Nant = length(positions)
    for antenna1 = 1:Nant, antenna2 = antenna1:Nant
        position1 = positions[antenna1]
        position2 = positions[antenna2]
        baseline = Baseline(baseline"ITRF",
                            position1.x - position2.x,
                            position1.y - position2.y,
                            position1.z - position2.z)
        push!(baselines, baseline)
    end
    baselines
end

