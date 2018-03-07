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

struct BlockVector{B} <: AbstractBlockVector
    path        :: String
    progressbar :: Bool
    distribute  :: Bool
    cached      :: Ref{Bool}
    size        :: Tuple{Int, Int}
    blocks      :: Vector{B}

    function BlockVector{B}(path, size, write=true;
                            progressbar=false, distribute=false, cached=false) where B
        if write
            isdir(path) || mkpath(path)
            save(joinpath(path, "METADATA.jld2"), "size", size)
        end
        blocks = Array{B}(size)
        output = new(path, progressbar, distribute, Ref(cached), size, blocks)
        if cached
            cache!(output)
        end
        output
    end
end

function BlockVector{B}(path; kwargs...) where B
    size = load(joinpath(path, "METADATA.jld2"), "size")
    BlockVector{B}(path, size, false; kwargs...)
end

function Base.show(io::IO, matrix::BlockVector)
    @printf(io, "BlockVector(%s, cached=%s)", matrix.path, matrix.cached[] ? "true" : "false")
end

indices(matrix::BlockVector) = [(idx, jdx) for jdx = 1:matrix.size[2] for idx = 1:matrix.size[1]]

function Base.getindex(matrix::BlockVector, idx, jdx)
    if matrix.cached[]
        return matrix.blocks[idx, jdx]
    else
        return read_from_disk(matrix, idx, jdx)
    end
end

function Base.setindex!(matrix::BlockVector, block, idx, jdx)
    if matrix.cached[]
        matrix.blocks[idx, jdx] = block
    else
        write_to_disk(matrix, block, idx, jdx)
    end
    block
end

function cache!(matrix::BlockVector)
    matrix.cached[] = true
    for jdx = 1:matrix.size[2], idx = 1:matrix.size[1]
        matrix.blocks[idx, jdx] = read_from_disk(matrix, idx, jdx)
    end
    matrix
end

function flush!(matrix::BlockVector)
    for jdx = 1:matrix.size[2], idx = 1:matrix.size[1]
        write_to_disk(matrix, matrix.blocks[idx, jdx], idx, jdx)
    end
    matrix.cached[] = false
    matrix
end

function read_from_disk(matrix::BlockVector{B}, idx, jdx) where B
    dirname    = @sprintf("%04d",      jdx)
    filename   = @sprintf("%04d.jld2", idx)
    objectname = "block"
    load(joinpath(matrix.path, dirname, filename), objectname) :: B
end

function write_to_disk(matrix::BlockVector{B}, block::B, idx, jdx) where B
    dirname    = @sprintf("%04d",      jdx)
    filename   = @sprintf("%04d.jld2", idx)
    objectname = "block"
    isdir(joinpath(matrix.path, dirname)) || mkpath(joinpath(matrix.path, dirname))
    save(joinpath(matrix.path, dirname, filename), objectname, block)
    block
end

