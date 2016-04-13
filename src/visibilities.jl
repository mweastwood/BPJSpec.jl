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

doc"""
    GriddedVisibilities

This type represents visibilities on a sidereal time grid.

# Fields

* `path` points to the directory where all of the visibilities are stored
* `Nbase` is the total number of baselines
* `Ntime` is the number of sidereal time grid points
* `frequencies` is the list of frequencies in units of Hz
* `origin` is the sidereal time of the first grid point
* `data` is a list of mmapped arrays where the visibilities are stored
* `weights` is a list of mmapped arrays where the weights are stored

Note that in each of `data` and `weights` there is one mmapped array for each
frequency channel, and each array has dimensions of `Nbase` by `Ntime`.

# Implementation

When gridding the visibilities we are going to be making a lot of small writes
to several arrays that may not fit into the system memory. This is why we mmap
the arrays onto the disk.

Also note that it seems like two arrays cannot be mmapped to the same file.
That is instead of being written one after another the two arrays are written
on top of each other. This is why `data` and `weights` are mmapped to two
separate files.
"""
immutable GriddedVisibilities
    path :: ASCIIString
    Nbase :: Int
    Ntime :: Int
    frequencies :: Vector{Float64}
    origin :: Float64
    data :: Vector{Matrix{Complex128}}
    weights :: Vector{Matrix{Float64}}
    function GriddedVisibilities(path, Nbase, Ntime, frequencies, origin, data, weights)
        0 ≤ origin < 1 || throw(ArgumentEttor("The sidereal time must be in the interval [0,1)"))
        new(path, Nbase, Ntime, frequencies, origin, data, weights)
    end
end

function GriddedVisibilities(path)
    local Nbase, Ntime, frequencies, origin
    # first read the METADATA file
    open(joinpath(path, "METADATA"), "r") do file
        Nbase = read(file, Int)
        Ntime = read(file, Int)
        len = read(file, Int)
        frequencies = read(file, Float64, len)
        origin = read(file, Float64)
    end
    # now mmap each frequency channel
    data = Matrix{Complex128}[]
    weights = Matrix{Float64}[]
    for frequency in frequencies
        filename = block_filename(frequency)
        open(joinpath(path, filename), "r+") do file
            push!(data, Mmap.mmap(file, Matrix{Complex128}, (Nbase, Ntime)))
        end
        open(joinpath(path, filename*".weights"), "r+") do file
            push!(weights, Mmap.mmap(file, Matrix{Float64}, (Nbase, Ntime)))
        end
    end
    GriddedVisibilities(path, Nbase, Ntime, frequencies, origin, data, weights)
end

function GriddedVisibilities(path, meta, Ntime)
    frequencies = meta.channels
    origin = sidereal_time(meta)
    # create the directory if it doesn't already exist
    isdir(path) || mkdir(path)
    # create the METADATA file to store Nbase, Ntime, and the list of frequency channels
    open(joinpath(path, "METADATA"), "w") do file
        write(file, Nbase(meta), Ntime, length(frequencies), frequencies, origin)
    end
    # create the files storying each frequency channel
    data = Matrix{Complex128}[]
    weights = Matrix{Float64}[]
    for channel = 1:length(frequencies)
        ν = frequencies[channel]
        filename = block_filename(ν)
        open(joinpath(path, filename), "w+") do file
            push!(data, Mmap.mmap(file, Matrix{Complex128}, (Nbase(meta), Ntime)))
        end
        open(joinpath(path, filename*".weights"), "w+") do file
            push!(weights, Mmap.mmap(file, Matrix{Float64}, (Nbase(meta), Ntime)))
        end
    end
    GriddedVisibilities(path, Nbase(meta), Ntime, frequencies, origin, data, weights)
end

Nfreq(visibilities::GriddedVisibilities) = length(visibilities.frequencies)

function getindex(visibilities::GriddedVisibilities, channel)
    data    = visibilities.data[channel]
    weights = visibilities.weights[channel]
    output  = similar(data)
    @inbounds for I in eachindex(data)
        output[I] = ifelse(weights[I] < eps(Float64), complex(0.0), data[I] / weights[I])
    end
    output
end

"""
    grid!(gridded_visibilities, meta, data)

Add the data to the sidereal time grid.
"""
function grid!(gridded_visibilities, meta, data)
    origin = gridded_visibilities.origin
    time = mod(sidereal_time(meta) - origin, 1)
    for idx1 = 1:Nfreq(meta)
        # The list of frequency channels we want to grid might not match
        # the list of frequency channels in the measurement set. Therefore
        # we need to skip over any frequencies that are not represented in
        # the gridded data.
        ν = meta.channels[idx1]
        idx2 = findfirst(gridded_visibilities.frequencies, ν)
        idx2 == 0 && continue
        # Now that we are sure that we want this frequency channel,
        # write this slice of the data onto the grid.
        grid = gridded_visibilities.data[idx2]
        weights = gridded_visibilities.weights[idx2]
        data_slice = slice(data.data, :, idx1)
        flag_slice = slice(data.flags, :, idx1)
        grid_onechannel!(grid, weights, data_slice, flag_slice, time)
    end
end

function grid_onechannel!(grid, weights, data, flags, time)
    Nbase, Ntime = size(grid)
    # Identify the two nearest time gridpoints
    # and calculate their corresponding weights
    Δt = 1/Ntime
    times = 0.0:Δt:(1.0-Δt)
    idx1 = searchsortedlast(times, time)
    idx2 = idx1 == Ntime? 1 : idx1+1
    weight1 = 1 - (time - times[idx1])/Δt
    weight2 = 1 - weight1 # note this handles the case where the second grid point has wrapped around
    # Write to the grid
    for α = 1:Nbase
        flags[α] && continue
        I = 0.5*(data[α].xx + data[α].yy)
        grid[α,idx1] += weight1 * I
        grid[α,idx2] += weight2 * I
        weights[α,idx1] += weight1
        weights[α,idx2] += weight2
    end
end

doc"""
    sidereal_time(meta)

Get the local sidereal time on the interval $0 \le t < 1$.
"""
function sidereal_time(meta)
    frame = TTCal.reference_frame(meta)
    zenith = Direction(dir"AZEL",0degrees,90degrees)
    zenith_app = measure(frame, zenith, dir"APP")
    time = longitude(zenith_app)
    mod2pi(time) / 2π
end

"""
    grid_visibilities(filename, ms::MeasurementSet)

Grid the data from the measurement set.
"""

doc"""
    grid_visibilities(filename, data, flags, frequencies, time)

For a sidereal time on the interval $0 \le t < 1$, grid and
add the data to the given file.

Here "grid" means to put the data onto a set of points that are
evenly space in sidereal time. If the correlator is continuously
running, it will spit out an integration every 30 seconds (for
example). However 30 seconds does not divide evenly into a sidereal
day so that these integrations will slowly start slipping with
respect to sidereal time as we begin to average more days together.
By gridding onto a sidereal time grid, we avoid this problem.
"""

