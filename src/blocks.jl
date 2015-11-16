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

abstract AbstractVectorBlock # block that composes a vector
abstract AbstractMatrixBlock # block that composes a matrix
typealias AbstractBlock Union{AbstractVectorBlock,AbstractMatrixBlock}

Base.size(A::AbstractBlock) = size(A.block)
Base.length(A::AbstractBlock) = length(A.block)
Base.svd(A::AbstractBlock) = svd(A.block,thin=true)

abstract AbstractBlockVector # a vector composed of blocks
abstract AbstractBlockDiagonalMatrix
abstract AbstractBlockMatrix # a matrix composed of blocks
typealias VectorOfBlocks Union{AbstractBlockVector,AbstractBlockDiagonalMatrix}

Nblocks(A::VectorOfBlocks) = length(A.blocks)

immutable MatrixBlock <: AbstractMatrixBlock
    block::Matrix{Complex128}
end

immutable BlockDiagonalMatrix <: AbstractBlockDiagonalMatrix
    blocks::Vector{MatrixBlock}
end

function *{T<:AbstractBlock}(lhs::AbstractMatrixBlock, rhs::T)
    meta   = metadata(rhs)
    result = lhs.block*rhs.block
    T(result,meta...)
end

function *{T<:VectorOfBlocks}(lhs::AbstractBlockDiagonalMatrix, rhs::T)
    Nblocks(lhs) == Nblocks(rhs) || error("Number of blocks must match.")
    T([lhs.blocks[i]*rhs.blocks[i] for i = 1:Nblocks(rhs)])
end

