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

struct BlockDiagonalMatrix <: BlockMatrix
    path :: String
    progressbar :: Bool
    distribute  :: Bool
    cached      :: Ref{Bool}
    mmax        :: Int
    blocks      :: Vector{Matrix{Complex128}}

    function BlockDiagonalMatrix(path, mmax, write=true;
                                 progressbar=false, distribute=false, cached=false)
        if write
            isdir(path) || mkpath(path)
            save(joinpath(path, "METADATA.jld2"), "mmax", mmax)
        end
        blocks = Matrix{Complex128}[]
        output = new(path, progressbar, distribute, Ref(cached), mmax, blocks)
        if cached
            cache!(output)
        end
        output
    end
end

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
    empty!(matrix.blocks)
    for m = 0:matrix.mmax
        push!(matrix.blocks, read_from_disk(matrix, m))
    end
    matrix
end

function flush!(matrix::BlockDiagonalMatrix)
    for m = 0:matrix.mmax
        write_to_disk(matrix, matrix.blocks[m+1], m)
    end
    empty!(matrix.blocks)
    matrix.cached[] = false
    matrix
end

function read_from_disk(matrix::BlockDiagonalMatrix, m)
    filename   = @sprintf("m=%04d.jld2", m)
    objectname = "block"
    load(joinpath(matrix.path, filename), objectname) :: Matrix{Complex128}
end

function write_to_disk(matrix::BlockDiagonalMatrix, block::Matrix{Complex128}, m)
    filename   = @sprintf("m=%04d.jld2", m)
    objectname = "block"
    save(joinpath(matrix.path, filename), objectname, block)
    block
end

