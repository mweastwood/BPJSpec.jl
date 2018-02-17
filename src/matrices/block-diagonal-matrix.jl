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

struct BlockDiagonalMatrix{B} <: BlockMatrix
    path        :: String
    progressbar :: Bool
    distribute  :: Bool
    cached      :: Ref{Bool}
    mmax        :: Int
    blocks      :: Vector{B}

    function BlockDiagonalMatrix{B}(path, mmax, write=true;
                                    progressbar=false, distribute=false, cached=false) where B
        if write
            isdir(path) || mkpath(path)
            save(joinpath(path, "METADATA.jld2"), "mmax", mmax)
        end
        blocks = Array{B}(mmax+1)
        output = new(path, progressbar, distribute, Ref(cached), mmax, blocks)
        if cached
            cache!(output)
        end
        output
    end
end

const DenseBlockDiagonalMatrix = BlockDiagonalMatrix{Matrix{Complex128}}

function BlockDiagonalMatrix(path; kwargs...)
    mmax = load(joinpath(path, "METADATA.jld2"), "mmax")
    BlockDiagonalMatrix(path, mmax, false; kwargs...)
end

Base.show(io::IO, matrix::BlockDiagonalMatrix) =
    print(io, "BlockDiagonalMatrix: ", matrix.path)

indices(matrix::BlockDiagonalMatrix) = collect(0:matrix.mmax)

function Base.getindex(matrix::BlockDiagonalMatrix, m)
    if matrix.cached[]
        return matrix.blocks[m+1]
    else
        return read_from_disk(matrix, m)
    end
end

function Base.setindex!(matrix::BlockDiagonalMatrix, block, m)
    if matrix.cached[]
        matrix.blocks[m+1] = block
    else
        write_to_disk(matrix, block, m)
    end
    block
end

function cache!(matrix::BlockDiagonalMatrix)
    matrix.cached[] = true
    for m = 0:matrix.mmax
        matrix.blocks[m+1] = read_from_disk(matrix, m)
    end
    matrix
end

function flush!(matrix::BlockDiagonalMatrix)
    for m = 0:matrix.mmax
        write_to_disk(matrix, matrix.blocks[m+1], m)
    end
    matrix.cached[] = false
    matrix
end

function read_from_disk(matrix::BlockDiagonalMatrix{B}, m) where B
    filename   = @sprintf("m=%04d.jld2", m)
    objectname = "block"
    load(joinpath(matrix.path, filename), objectname) :: B
end

function write_to_disk(matrix::BlockDiagonalMatrix{B}, block::B, m) where B
    filename   = @sprintf("m=%04d.jld2", m)
    objectname = "block"
    save(joinpath(matrix.path, filename), objectname, block)
    block
end

