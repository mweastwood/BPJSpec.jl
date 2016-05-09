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

function MModes(path, mmax, frequencies)
    # create the directory if it doesn't already exist
    isdir(path) || mkdir(path)
    # create the METADATA file to store mmax and the list of frequency channels
    open(joinpath(path, "METADATA"), "w") do file
        write(file, mmax, length(frequencies), frequencies)
    end
    # create the files storing each block
    blocks = Vector{Complex128}[]
    for frequency in frequencies, m = 0:mmax
        filename = block_filename(m, frequency)
        open(joinpath(path, filename), "w+") do file
            # the only sensible default here is to set the length of each block
            # to zero initially because we don't yet know how long each block
            # actually needs to be (after compression the length of each block
            # can be anything)
            len = 0
            write(file, len)
            # we can't mmap a zero-length array so just push an empty array
            push!(blocks, zeros(Complex128, len))
        end
    end
    MModes(path, mmax, frequencies, blocks)
end

function MModes(path)
    local mmax, frequencies
    # first read the METADATA file
    open(joinpath(path, "METADATA"), "r") do file
        mmax = read(file, Int)
        len  = read(file, Int)
        frequencies = read(file, Float64, len)
    end
    # now mmap each block
    blocks = Vector{Complex128}[]
    for frequency in frequencies, m = 0:mmax
        filename = block_filename(m, frequency)
        open(joinpath(path, filename), "r+") do file
            len = read(file, Int)
            push!(blocks, Mmap.mmap(file, Vector{Complex128}, len))
        end
    end
    MModes(path, mmax, frequencies, blocks)
end

function MModes(path, transfermatrix::TransferMatrix, alm::Alm)
    # assume the transfer matrix is single frequency for now
    mmodes = MModes(path, transfermatrix.mmax, transfermatrix.frequencies)
    for m = 0:transfermatrix.mmax
        @show m
        @time Bm = transfermatrix[m,1]
        @time am = getblock(alm, m)
        @time vm = Bm*am
        @time ν = mmodes.frequencies[1]
        @time mmodes[m,1] = vm
    end
end

doc"""
    MModes(path, visibilities::GriddedVisibilities, mmax)

Calculate the $m$-modes from the given visibilities. The visibilities
should be provided as a matrix where the first dimension indicates
the baseline and the second dimension indicates the sidereal time.
The FFT will be performed over the second dimension. The
visibilities should span a full sidereal day.
"""
function MModes(path, visibilities::GriddedVisibilities, mmax)
    mmodes = MModes(path, mmax, visibilities.frequencies)
    # create the files for storing each block
    for channel = 1:Nfreq(mmodes), m = 0:mmax
        mmodes[m,channel] = zeros(Complex128, two(m)*visibilities.Nbase)
    end
    # Fourier transform the visibilities with respect to sidereal time and pack
    # them into a vector
    for channel = 1:Nfreq(mmodes)
        info("Reading visibilities")
        grid = visibilities[channel]
        info("Performing Fourier transform")
        transformed_visibilities = do_fourier_transform(grid)
        info("Writing m-modes")
        pack_mmodes!(mmodes.blocks, transformed_visibilities, mmax, channel)
    end
    mmodes
end

function do_fourier_transform(matrix)
    Nbase, Ntime = size(matrix)
    planned_fft = plan_fft(matrix, 2)
    planned_fft * matrix / Ntime
end

function pack_mmodes!(blocks, transformed_visibilities, mmax, channel)
    Nbase, Ntime = size(transformed_visibilities)
    for m = 0:mmax
        idx = block_index(mmax, m, channel)
        block = blocks[idx]
        if m == 0
            for α = 1:Nbase
                block[α] = transformed_visibilities[α,1]
            end
        else
            for α = 1:Nbase
                α1 = 2α-1 # positive m
                α2 = 2α-0 # negative m
                block[α1] =      transformed_visibilities[α,m+1]
                block[α2] = conj(transformed_visibilities[α,Ntime+1-m])
            end
        end
    end
end

Nfreq(mmodes::MModes) = length(mmodes.frequencies)

function getindex(mmodes::MModes, m, channel)
    idx = block_index(mmodes.mmax, m, channel)
    copy(mmodes.blocks[idx])
end

function setindex!(mmodes::MModes, block, m, channel)
    idx = block_index(mmodes.mmax, m, channel)
    frequency = mmodes.frequencies[channel]
    filename = block_filename(m, frequency)
    open(joinpath(mmodes.path, filename), "w+") do file
        finalize(mmodes.blocks[idx]) # delete the previous mmapped array first
        write(file, length(block))
        mmodes.blocks[idx] = Mmap.mmap(file, Vector{Complex128}, length(block))
    end
    mmodes.blocks[idx][:] = block
end

