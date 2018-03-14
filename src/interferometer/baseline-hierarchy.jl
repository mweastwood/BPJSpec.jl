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
    Nfreq     :: Int
end

function Base.show(io::IO, hierarchy::Hierarchy)
    total_space = 0.0
    total_baselines = 0
    @printf(io, "|--------------+-----------+-------------|\n")
    @printf(io, "| lmin to lmax | baselines |  disk space |\n")
    @printf(io, "|--------------+-----------+-------------|\n")
    for idx = 1:length(hierarchy.divisions)-1
        lmin = hierarchy.divisions[idx]
        lmax = hierarchy.divisions[idx+1]
        Nbase = length(hierarchy.baselines[idx])
        space = Nbase*lmax*lmax*128/(8*1024^3)*hierarchy.Nfreq
        @printf(io, "| %4d to %4d |     %5d | %8.3f GB |\n", lmin, lmax, Nbase, space)
        total_space += space
        total_baselines += Nbase
    end
    @printf(io, "|--------------+-----------+-------------|\n")
    @printf(io, "|        total |     %5d | %8.3f GB |\n", total_baselines, total_space)
    @printf(io, "|--------------+-----------+-------------|\n")
end

function compute_baseline_hierarchy(metadata::Metadata, cutoff)
    lmax = maximum_multipole_moment(metadata)
    divisions = identify_divisions(lmax, cutoff)
    divisions = remove_empty_intervals(divisions)
    divisions = establish_minimum_lmax(divisions)
    baselines = categorize_baselines(lmax, divisions)
    Hierarchy(divisions, baselines, length(metadata.frequencies))
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
    lmax_range = 0:max(maximum(lmax), cutoff)+1
    histogram = zeros(length(lmax_range))
    for l in lmax
        histogram[l+1] += 1
    end
    histogram = histogram[1:cutoff+1]
    cumulative_histogram = cumsum(histogram)
    depth = 4
    find_crossover_points(cumulative_histogram, depth)
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

function remove_empty_intervals(divisions)
    unique(divisions)
end

function establish_minimum_lmax(divisions)
    # There are accuracy issues associated with computing spherical harmonics with a small lmax.
    # This comes about because the resolution of the map scales with lmax and so the beam is poorly
    # sampled when lmax is small. In the future we should apply some padding to ther spherical
    # harmonic transform such that this is not so much of an issue any more, but for now we'll just
    # increase lmax for the shortest baselines.
    minimum_lmax = 32
    idx = searchsortedlast(divisions, minimum_lmax)
    [[0, 32]; divisions[idx+1:end]]
end

"""
Categorize each baseline based on the divisions computed above.
"""
function categorize_baselines(lmax, divisions)
    [find(divisions[idx] .≤ lmax .< divisions[idx+1]) for idx = 1:length(divisions)-1]
end

function Nbase(hierarchy::Hierarchy, l)
    output = 0
    for idx = 1:length(hierarchy.divisions)-1
        lmax = hierarchy.divisions[idx+1]
        l > lmax && continue
        output += length(hierarchy.baselines[idx])
    end
    output
end

