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

function GriddedVisibilities(path, Nbase, Ntime, frequencies, origin)
    # create the directory if it doesn't already exist
    isdir(path) || mkdir(path)
    # create the METADATA file to store Nbase, Ntime, and the list of frequency channels
    open(joinpath(path, "METADATA"), "w") do file
        write(file, Nbase, Ntime, length(frequencies), frequencies, origin)
    end
    # create the files storying each frequency channel
    data = Matrix{Complex128}[]
    weights = Matrix{Float64}[]
    for frequency in frequencies
        filename = block_filename(frequency)
        open(joinpath(path, filename), "w+") do file
            push!(data, Mmap.mmap(file, Matrix{Complex128}, (Nbase, Ntime)))
        end
        open(joinpath(path, filename*".weights"), "w+") do file
            push!(weights, Mmap.mmap(file, Matrix{Float64}, (Nbase, Ntime)))
        end
    end
    GriddedVisibilities(path, Nbase, Ntime, frequencies, origin, data, weights)
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

function GriddedVisibilities(path, meta::Metadata, Ntime)
    frequencies = meta.channels
    origin = sidereal_time(meta)
    GriddedVisibilities(path, Nbase(meta), Ntime, frequencies, origin)
end

"""
    GriddedVisibilities(path, meta::Metadata, mmodes::MModes, Ntime)

Calculate the visibilities from the given m-modes.
"""
function GriddedVisibilities(path, meta::Metadata, mmodes::MModes, Ntime)
    origin = sidereal_time(meta)
    visibilities = GriddedVisibilities(path, Nbase(meta), Ntime, mmodes.frequencies, origin)
    for idx = 1:Nfreq(visibilities)
        block = zeros(Complex128, Nbase(meta), Ntime)
        # m = 0
        v = mmodes[0,idx]
        for α = 1:Nbase(meta)
            block[α,1] = v[α]
        end
        # m > 0
        for m = 1:mmodes.mmax
            v = mmodes[m,idx]
            for α = 1:Nbase(meta)
                α1 = 2α-1 # positive m
                α2 = 2α-0 # negative m
                block[α,m+1]       =      v[α1]
                block[α,Ntime+1-m] = conj(v[α2])
            end
        end
        info("Starting FFT...")
        @time result = ifft(block,2)
        info("Doing the multiply...")
        @time result = result*Ntime
        info("Finished FFT...")
        visibilities.data[idx][:] = result
        visibilities.weights[idx][:] = 1
    end
    visibilities
end

Nfreq(visibilities::GriddedVisibilities) = length(visibilities.frequencies)

function getindex(visibilities::GriddedVisibilities, channel)
    #data    = visibilities.data[channel]
    #weights = visibilities.weights[channel]
    local data, weights
    filename = block_filename(visibilities.frequencies[channel])
    open(joinpath(visibilities.path, filename), "r") do file
        data = read(file, Complex128, (visibilities.Nbase, visibilities.Ntime))
    end
    open(joinpath(visibilities.path, filename*".weights"), "r") do file
        weights = read(file, Float64, (visibilities.Nbase, visibilities.Ntime))
    end

    output  = similar(data)
    @inbounds for j = 1:visibilities.Ntime, i = 1:visibilities.Nbase
        d = data[i,j]
        w = weights[i,j]
        output[i,j] = ifelse(w < eps(Float64), complex(0.0), d/w)
    end
    output
end

function setindex!(visibilities::GriddedVisibilities, data, channel)
    weights = ones(Float64, size(data))
    filename = block_filename(visibilities.frequencies[channel])
    open(joinpath(visibilities.path, filename), "w") do file
        write(file, data)
    end
    open(joinpath(visibilities.path, filename*".weights"), "w") do file
        write(file, weights)
    end
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
    zenith = Direction(dir"AZEL", 0degrees, 90degrees)
    zenith_app = measure(frame, zenith, dir"APP")
    time = longitude(zenith_app)
    mod2pi(time) / 2π
end

