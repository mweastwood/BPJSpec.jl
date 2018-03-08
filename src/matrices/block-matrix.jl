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

doc"""
    struct BlockMatrix{B, M, S, N}

This type represents a (potentially enormous) block-diagonal matrix. This type is designed to be
general enough to handle large matrices that fit in memory as well as enormous matrices that do not
fit in memory.

# Type Parameters

* `B` specifies the type of the blocks that compose the matrix
* `N` specifies the number of indices used to index into the matrix blocks
* `M` specifies the type of metadata that is associated with the matrix
* `S` specifies the type of storage used for the matrix

# Fields


"""
struct BlockMatrix{B, N, M, S}
    metadata :: M
    storage  :: S
    cache    :: Cache{B, N}
end

function BlockMatrix{B, N}(storage::Mechanism, metadata::MatrixMetadata, write=true) where {B, N}
    write && write_metadata(storage, metadata)
    BlockMatrix(metadata, storage, Cache{B}(number_of_blocks(metadata)))
end

function BlockMatrix{B, N}(storage::NoFile, metadata::MatrixMetadata) where {B, N}
    cache = set!(Cache{B}(number_of_blocks(metadata)))
    BlockMatrix(metadata, storage, cache)
end

function BlockMatrix{B, N}(metadata::MatrixMetadata) where {B, N}
    BlockMatrix{B, N}(NoFile(), metadata)
end

function BlockMatrix{B, N}(path::String) where {B, N}
    storage, metadata = read_metadata(path)
    BlockMatrix{B, N}(storage, metadata, false)
end

function Base.show(io::IO, matrix::BlockMatrix)
    @printf(io, "BlockMatrix(%s, metadata=%s, cached=%s)",
            matrix.storage, matrix.metadata, used(matrix.cache) ? "true" : "false")
end

indices(matrix::BlockMatrix{B, 1}) where {B} = 1:number_of_blocks(matrix.metadata)

function indices(matrix::BlockMatrix{B, 2}) where {B}
    N, M = number_of_blocks(matrix.metadata)
    ((idx, jdx) for jdx = 1:M for idx = 1:N)
end

function Base.getindex(matrix::BlockMatrix{B}, idx...)::B where B
    if used(matrix.cache)
        return matrix.cache[idx...]
    else
        return matrix.storage[idx...]
    end
end

function Base.setindex!(matrix::BlockMatrix{B}, block::B, idx...) where B
    if used(matrix.cache)
        matrix.cache[idx...] = block
    else
        matrix.storage[idx...] = block
    end
    block
end

function cache!(matrix::BlockMatrix)
    set!(matrix.cache)
    for idx in indices(matrix)
        matrix.cache[idx] = matrix.storage[idx]
    end
    matrix
end

function flush!(matrix::BlockMatrix)
    unset!(matrix.cache)
    for idx in indices(matrix)
        matrix.storage[idx] = matrix.cache[idx]
    end
    matrix
end

distribute_write(matrix::BlockMatrix) = distribute_write(matrix.storage) && !used(matrix.cache)
distribute_read(matrix::BlockMatrix) = distribute_read(matrix.storage) && !used(matrix.cache)

