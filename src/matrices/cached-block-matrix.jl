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

struct CachedBlockMatrix{M, B} <: BlockMatrix
    matrix :: M
    blocks :: B
end

function cache(matrix::BlockMatrix)
    blocks = [matrix[idx...] for idx in Iterators.product(indices(matrix)...)]
    CachedBlockMatrix(matrix, blocks)
end

Base.indices(matrix::CachedBlockMatrix) = indices(matrix.matrix)
normalize_indices(::BlockMatrix, idx) = idx
normalize_indices(::BlockMatrix, idx, jdx) = (idx, jdx)
@inline normalize_indices(::BlockMatrix, idx...) = idx

@inline function Base.getindex(matrix::CachedBlockMatrix, idx...)
    idx = normalize_indices(matrix.matrix, idx...)
    matrix.blocks[idx...]
end

@inline function Base.setindex(matrix::CachedBlockMatrix, block, idx...)
    idx = normalize_indices(matrix.matrix, idx...)
    matrix.blocks[idx...] = block
end

function flush(matrix::CachedBlockMatrix)
    for idx in Iterators.product(indices(matrix))
        matrix.matrix[idx...] = matrix[idx]
    end
end

