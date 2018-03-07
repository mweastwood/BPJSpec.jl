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

struct BlockMatrix{B, M, S} <: AbstractBlockMatrix
    metadata :: M
    storage  :: S
    cache    :: Cache{B, 2}
end

function BlockMatrix{B}(storage::Mechanism, metadata::MatrixMetadata, write=true) where B
    write && write_metadata(storage, metadata)
    BlockMatrix(metadata, storage, Cache{B}(number_of_blocks(metadata)...))
end

function BlockMatrix{B}(storage::NoFile, metadata::MatrixMetadata) where B
    cache = set!(Cache{B}(number_of_blocks(metadata)...))
    BlockMatrix(metadata, storage, cache)
end

function BlockMatrix{B}(metadata::MatrixMetadata) where B
    BlockMatrix{B}(NoFile(), metadata)
end

function BlockMatrix{B}(path::String) where B
    storage, metadata = read_metadata(path)
    BlockMatrix{B}(storage, metadata, false)
end

function Base.show(io::IO, matrix::BlockMatrix)
    @printf(io, "BlockMatrix(%s, metadata=%s, cached=%s)",
            matrix.path, matrix.metadata, used(matrix.cache) ? "true" : "false")
end

function indices(matrix::BlockMatrix)
    N, M = number_of_blocks(matrix.metadata)
    ((idx, jdx) for jdx = 1:M for idx = 1:N)
end

function Base.getindex(matrix::BlockMatrix{B}, idx, jdx)::B where B
    if used(matrix.cache)
        return matrix.cache[idx, jdx]
    else
        return matrix.storage[idx, jdx]
    end
end

function Base.setindex!(matrix::BlockMatrix{B}, block::B, idx, jdx) where B
    if used(matrix.cache)
        matrix.cache[idx, jdx] = block
    else
        matrix.storage[idx, jdx] = block
    end
    block
end

function cache!(matrix::BlockMatrix)
    set!(matrix.cache)
    N, M = number_of_blocks(matrix.metadata)
    for jdx = 1:M, idx = 1:N
        matrix.cache[idx, jdx] = matrix.storage[idx, jdx]
    end
    matrix
end

function flush!(matrix::BlockMatrix)
    unset!(matrix.cache)
    N, M = number_of_blocks(matrix.metadata)
    for jdx = 1:M, idx = 1:N
        matrix.storage[idx, jdx] = matrix.cache[idx, jdx]
    end
    matrix
end

