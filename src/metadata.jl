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

function baseline_hierarchy(metadata::Metadata)
    νmax = maximum(metadata.frequencies)
    λmin = u"c" / νmax
    lmax_float = Float64.(2π .* norm.(metadata.baselines) ./ λmin)
    # Add some extra slack because the sensitivity of a baseline of a given length extends a little
    # bit past the edge of the estimate given on the previous line.
    lmax = ceil.(Int, lmax_float .* 1.1)

    # Compute the category boundaries
    lmax_range = 0:maximum(lmax)+1
    histogram = zeros(length(lmax_range))
    for l in lmax
        histogram[l+1] += 1
    end
    cumulative_histogram = cumsum(histogram)
    divisions = find_crossover_points(cumulative_histogram, 4)

    # Separate the baselines into the categories
    categories = zeros(Int, length(metadata.baselines))
    for idx = 1:length(divisions)-1
        categories[divisions[idx] .≤ lmax .< divisions[idx+1]] = idx
    end

    lmax, divisions, categories, cumulative_histogram
end

function find_crossover_points(cumulative_histogram, depth)
    lmin = 0
    lmax = length(cumulative_histogram) - 1
    divisions = [lmin, lmax]
    find_crossover_points!(divisions, cumulative_histogram, lmin, lmax, depth)
    sort!(divisions)
end

function find_crossover_points!(divisions, cumulative_histogram, lmin, lmax, depth)
    if depth > 0
        l = find_crossover_point(cumulative_histogram, lmin, lmax)
        push!(divisions, l)
        find_crossover_points!(divisions, cumulative_histogram, lmin, l+0, depth-1)
        find_crossover_points!(divisions, cumulative_histogram, l+1, lmax, depth-1)
    end
end

function find_crossover_point(cumulative_histogram, lmin, lmax)
    for l = lmin:lmax
        space_below = (cumulative_histogram[l+1]-cumulative_histogram[lmin+1]) * (l+1)^2
        space_above = (cumulative_histogram[lmax+1]-cumulative_histogram[l+1]) * (lmax+1)^2
        if space_below > space_above
            return l
        end
    end
    lmax
end

# compatibility with TTCal

function from_ttcal(ttcal_metadata)
    frequencies  = ttcal_metadata.channels * u"Hz"
    position     = ttcal_position(ttcal_metadata)
    baselines    = ttcal_baselines(ttcal_metadata)
    phase_center = ttcal_metadata.phase_center
    Metadata(frequencies, position, baselines, phase_center)
end

function ttcal_position(ttcal_metadata)
    antenna_positions = getfield.(ttcal_metadata.antennas, :position)
    mean(antenna_positions)
end

function ttcal_baselines(ttcal_metadata)
    antenna_positions = getfield.(ttcal_metadata.antennas, :position)
    baselines = Baseline[]
    for α = 1:length(ttcal_metadata.baselines)
        antenna1 = antenna_positions[ttcal_metadata.baselines[α].antenna1]
        antenna2 = antenna_positions[ttcal_metadata.baselines[α].antenna2]
        baseline = Baseline(baseline"ITRF",
                            antenna1.x - antenna2.x,
                            antenna1.y - antenna2.y,
                            antenna1.z - antenna2.z)
        push!(baselines, baseline)
    end
    baselines
end

