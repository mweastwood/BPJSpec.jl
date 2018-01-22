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

struct Hierarchy
    divisions :: Vector{Int}
    baselines :: Vector{Vector{Int}}
end

function Base.show(io::IO, hierarchy::Hierarchy)
    total_space = 0.0
    total_baselines = 0
    @printf("|--------------+-----------+-------------|\n")
    @printf("| lmin to lmax | baselines |  disk space |\n")
    @printf("|--------------+-----------+-------------|\n")
    for idx = 1:length(hierarchy.divisions)-1
        lmin = hierarchy.divisions[idx]
        lmax = hierarchy.divisions[idx+1]
        Nbase = length(hierarchy.baselines[idx])
        space = Nbase*lmax*lmax*128/(8*1024^3)
        @printf("| %4d to %4d |     %5d | %8.3f GB |\n", lmin, lmax, Nbase, space)
        total_space += space
        total_baselines += Nbase
    end
    @printf("|--------------+-----------+-------------|\n")
    @printf("|        total |     %5d | %8.3f GB |\n", total_baselines, total_space)
    @printf("|--------------+-----------+-------------|\n")
end

function compute_baseline_hierarchy(metadata::Metadata, cutoff)
    lmax = maximum_multipole_moment(metadata)
    divisions = identify_divisions(lmax, cutoff)
    baselines = categorize_baselines(lmax, divisions)
    Hierarchy(divisions, baselines)
end

"""
Compute the maximum multipole moment l that each baseline is sensitive to.
"""
function maximum_multipole_moment(metadata)
    νmax = maximum(metadata.frequencies)
    λmin = u"c" / νmax
    lmax_float = Float64.(2π .* norm.(metadata.baselines) ./ λmin)
    # Add some extra slack because the sensitivity of a baseline of a given length extends a little
    # bit past the edge of the estimate given on the previous line.
    ceil.(Int, lmax_float .* 1.1)
end

"""
Separate the baselines based on their length. We will try to do this by recursively splitting the
baselines into two groups where each group requires roughly the same amount of space (in the
transfer matrix).
"""
function identify_divisions(lmax, cutoff)
    lmax_range = 0:maximum(lmax)+1
    histogram = zeros(length(lmax_range))
    for l in lmax
        histogram[l+1] += 1
    end
    histogram = histogram[1:cutoff+1]
    cumulative_histogram = cumsum(histogram)
    find_crossover_points(cumulative_histogram, 4)
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
        # storage space required to store the baselines shorter than l in the transfer matrix
        space_below = (cumulative_histogram[l+1]-cumulative_histogram[lmin+1]) * (l+1)^2
        # storage space required to store the baselines longer than l in the transfer matrix
        space_above = (cumulative_histogram[lmax+1]-cumulative_histogram[l+1]) * (lmax+1)^2
        if space_below > space_above
            return l
        end
    end
    lmax
end

"""
Categorize each baseline based on the divisions computed above.
"""
function categorize_baselines(lmax, divisions)
    [find(divisions[idx] .≤ lmax .< divisions[idx+1]) for idx = 1:length(divisions)-1]
end

