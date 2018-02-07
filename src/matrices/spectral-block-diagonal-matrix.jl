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

struct SpectralBlockDiagonalMatrix <: BlockMatrix
    path        :: String
    progressbar :: Bool
    distribute  :: Bool
    cached      :: Ref{Bool}
    mmax        :: Int
    frequencies :: Vector{typeof(1.0*u"Hz")}
    blocks      :: Matrix{Matrix{Complex128}}

    function SpectralBlockDiagonalMatrix(path, mmax, frequencies, write=true;
                                         progressbar=false, distribute=false, cached=false)
        if write
            isdir(path) || mkpath(path)
            save(joinpath(path, "METADATA.jld2"), "mmax", mmax, "frequencies", frequencies)
        end
        blocks = Array{Matrix{Complex128}}(mmax+1, length(frequencies))
        output = new(path, progressbar, distribute, Ref(cached), mmax, frequencies, blocks)
        if cached
            cache!(output)
        end
        output
    end
end

function SpectralBlockDiagonalMatrix(path; kwargs...)
    mmax, frequencies = load(joinpath(path, "METADATA.jld2"), "mmax", "frequencies")
    SpectralBlockDiagonalMatrix(path, mmax, frequencies, false; kwargs...)
end

Base.show(io::IO, matrix::SpectralBlockDiagonalMatrix) =
    print(io, "SpectralBlockDiagonalMatrix: ", matrix.path)

indices(matrix::SpectralBlockDiagonalMatrix) =
    [(m, β) for β = 1:length(matrix.frequencies) for m = 0:matrix.mmax]

function Base.getindex(matrix::SpectralBlockDiagonalMatrix, m, β::Integer)
    if matrix.cached[]
        return matrix.blocks[m+1, β]
    else
        return read_from_disk(matrix, m, β)
    end
end

function Base.setindex!(matrix::SpectralBlockDiagonalMatrix, block, m, β)
    if matrix.cached[]
        matrix.blocks[m+1, β] = block
    else
        write_to_disk(matrix, block, m, β)
    end
    block
end

function Base.getindex(matrix::SpectralBlockDiagonalMatrix, m, β::AbstractVector)
    blocks = [matrix[m, β′] for β′ in β]
    X = sum(size.(blocks, 1))
    Y = sum(size.(blocks, 2))
    output = zeros(Complex128, X, Y)
    x = y = 1
    for block in blocks
        output[x:x+size(block, 1)-1, y:y+size(block, 2)-1] = block
        x += size(block, 1)
        y += size(block, 2)
    end
    output
end

Base.getindex(matrix::SpectralBlockDiagonalMatrix, m) = matrix[m, 1:length(matrix.frequencies)]

"""
Stack the given frequency blocks on top of each other (as opposed to diagonally). This is used when
averaging frequency channels together.
"""
function stack(matrix::SpectralBlockDiagonalMatrix, m, β::AbstractVector)
    blocks = [matrix[m, β′] for β′ in β]
    X = sum(size.(blocks, 1))
    Y = size(blocks[1], 2)
    output = zeros(Complex128, X, Y)
    x = 1
    for block in blocks
        output[x:x+size(block, 1)-1, :] = block
        x += size(block, 1)
    end
    output
end

function cache!(matrix::SpectralBlockDiagonalMatrix)
    matrix.cached[] = true
    for β = 1:length(matrix.frequencies), m = 0:matrix.mmax
        matrix.blocks[m+1, β] = read_from_disk(matrix, m, β)
    end
    matrix
end

function flush!(matrix::SpectralBlockDiagonalMatrix)
    for β = 1:length(matrix.frequencies), m = 0:matrix.mmax
        write_to_disk(matrix, matrix.blocks[m+1, β], m, β)
    end
    matrix.cached[] = false
    matrix
end

function read_from_disk(matrix::SpectralBlockDiagonalMatrix, m, β)
    ν = matrix.frequencies[β]
    dirname    = @sprintf("%.3fMHz", ustrip(uconvert(u"MHz", ν)))
    filename   = @sprintf("%04d.jld2", m)
    objectname = "block"
    load(joinpath(matrix.path, dirname, filename), objectname) :: Matrix{Complex128}
end

function write_to_disk(matrix::SpectralBlockDiagonalMatrix, block::Matrix{Complex128}, m, β)
    ν = matrix.frequencies[β]
    dirname    = @sprintf("%.3fMHz", ustrip(uconvert(u"MHz", ν)))
    filename   = @sprintf("%04d.jld2", m)
    objectname = "block"
    isdir(joinpath(matrix.path, dirname)) || mkpath(joinpath(matrix.path, dirname))
    save(joinpath(matrix.path, dirname, filename), objectname, block)
    block
end

