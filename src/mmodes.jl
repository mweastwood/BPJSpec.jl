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
    MModes

This type represents the $m$-modes measured by an interferometer.

$m$-modes are the Fourier transform of a visibility with
respect to sidereal time. These are related to the spherical
harmonic coefficients of the sky brightness through a matrix equation.

\[
    v = Ba,
\]

where $v$ is the vector of $m$-modes, $B$ is the transfer
matrix, and $a$ is the vector of spherical harmonic coefficients.

# Fields

* `path` points to the directory that contains the blocks of $m$-modes
* `mmax` is the maximum value of the azimuthal quantum number $m$
* `frequencies` is the list of frequencies in units of Hz
* `blocks` is a list of mmapped vectors where the actual $m$-modes are stored

Note that there is one mmapped vector for each frequency and each value of $m$.
"""
immutable MModes
    path :: ASCIIString
    mmax :: Int
    frequencies :: Vector{Float64}
    blocks :: Vector{Vector{Complex128}}
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

doc"""
    MModes(path, visibilities, mmax)

Calculate the $m$-modes from the given visibilities. The visibilities
should be provided as a matrix where the first dimension indicates
the baseline and the second dimension indicates the sidereal time.
The FFT will be performed over the second dimension. The
visibilities should span a full sidereal day.
"""
function MModes(path, visibilities, mmax)
    frequencies = visibilities.frequencies
    # create the directory if it doesn't already exist
    isdir(path) || mkdir(path)
    # create the METADATA file to store mmax and the list of frequency channels
    open(joinpath(path, "METADATA"), "w") do file
        write(file, mmax, length(frequencies), frequencies)
    end
    # create the files for storing each block
    blocks = Vector{Complex128}[]
    for channel = 1:length(frequencies), m = 0:mmax
        ν = frequencies[channel]
        filename = block_filename(m, ν)
        open(joinpath(path, filename), "w+") do file
            len = two(m)*visibilities.Nbase
            write(file, len)
            push!(blocks, Mmap.mmap(file, Vector{Complex128}, len))
        end
    end
    # Fourier transform the visibilities with respect to sidereal time and pack
    # them into a vector
    for channel = 1:length(frequencies)
        transformed_visibilities = do_fourier_transform(visibilities[channel])
        pack_mmodes!(blocks, transformed_visibilities, mmax, channel)
    end
    MModes(path, mmax, frequencies, blocks)
end

function do_fourier_transform(matrix)
    Nbase, Ntime = size(matrix)
    fft(matrix, 2) / Ntime
end

function pack_mmodes!(blocks, transformed_visibilities, mmax, channel)
    Nbase, Ntime = size(transformed_visibilities)
    for m = 0:mmax
        idx = block_index(mmax, m, channel)
        block = blocks[idx]
        if m == 0
            for α = 1:Nbase
                block[α] = transformed_visibilities[α,m+1]
            end
        else
            for α = 1:Nbase
                α1 = α         # positive m
                α2 = α + Nbase # negative m
                block[α1] =      transformed_visibilities[α,m+1]
                block[α2] = conj(transformed_visibilities[α,Ntime+1-m])
            end
        end
    end
end

Nfreq(mmodes::MModes) = length(mmodes.ν)

function getindex(mmodes::MModes, m, channel)
    idx = block_index(mmodes.mmax, m, channel)
    copy(mmodes.blocks[idx])
end

#=
"""
    visibilities(v::MModes)

Calculate the visibilities from the given m-modes. The visibilities
will be returned as a two dimensional array where the first dimension
indicates the baseline and the second dimension indicates the sidereal
time, which evenly spans a full sidereal day.
"""
function visibilities(v::MModes)
    is_single_frequency(v) || error("Expected single-frequency m-modes.")
    Nbase = size(v[1].block,1)
    Ntime = 2mmax(v)+1
    M = zeros(Complex128,Nbase,Ntime)
    for α = 1:Nbase
        M[α,1] = v[1][α]
    end
    for m = 1:mmax(v), α = 1:Nbase
        α1 = α         # positive m
        α2 = α + Nbase # negative m
        M[α,m+1]       =      v[m+1][α1]
        M[α,Ntime+1-m] = conj(v[m+1][α2])
    end
    ifft(M,2)*Ntime
end
=#

#=
function call(::Type{MModes}, vector::Vector{Complex128}, m::Int, ν::AbstractVector)
    # When foreground-filtering the m-modes, it is often convenient to call
    # `full(::MModes)` to turn the m-modes into a standard vector. This function
    # is the inverse of `full`.
    Nfreq = length(ν)
    N = div(length(vector),Nfreq)
    blocks = [VectorBlock(vector[(β-1)*N+1:β*N]) for β = 1:Nfreq]
    meta = MModesMeta(m:m,collect(ν))
    MModes(blocks,meta)
end
=#

