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

abstract AbstractBlockVector # a vector composed of blocks
abstract AbstractBlockDiagonalMatrix
typealias ListOfBlocks Union{AbstractBlockVector,AbstractBlockDiagonalMatrix}

metadata(::AbstractBlock) = ()

getindex(A::AbstractVectorBlock,i) = A.block[i]
setindex!(A::AbstractVectorBlock,x,i) = A.block[i] = x
getindex(A::AbstractMatrixBlock,i,j) = A.block[i,j]
setindex!(A::AbstractMatrixBlock,x,i,j) = A.block[i,j] = x

Base.size(A::AbstractBlock) = size(A.block)
Base.length(A::AbstractBlock) = length(A.block)

function Base.length(A::AbstractBlockVector)
    sz = 0
    for block in A.blocks
        sz += length(block)
    end
    sz
end

function Base.size(A::AbstractBlockDiagonalMatrix)
    sz1 = sz2 = 0
    for block in A.blocks
        sz = size(block)
        sz1 += sz[1]
        sz2 += sz[2]
    end
    sz1,sz2
end

Nblocks(A::ListOfBlocks) = length(A.blocks)

==(lhs::AbstractBlock, rhs::AbstractBlock) = metadata(lhs) == metadata(rhs) && lhs.block == rhs.block
==(lhs::ListOfBlocks, rhs::ListOfBlocks) = lhs.blocks == rhs.blocks

function Base.full(A::AbstractBlockVector)
    B = zeros(Complex128,length(A))
    idx = 1
    for block in A.blocks
        l = length(block)
        B[idx:idx+l-1] = block.block
        idx += l
    end
    B
end

function Base.full(A::AbstractBlockDiagonalMatrix)
    B = zeros(Complex128,size(A))
    idx1 = idx2 = 1
    for block in A.blocks
        sz = size(block)
        B[idx1:idx1+sz[1]-1,idx2:idx2+sz[2]-1] = block.block
        idx1 += sz[1]
        idx2 += sz[2]
    end
    B
end

function *{T<:AbstractBlock}(lhs::AbstractMatrixBlock, rhs::T)
    meta   = metadata(rhs)
    result = lhs.block*rhs.block
    T(result,meta...)
end

function *{T<:ListOfBlocks}(lhs::AbstractBlockDiagonalMatrix, rhs::T)
    Nblocks(lhs) == Nblocks(rhs) || error("Number of blocks must match.")
    T([lhs.blocks[i]*rhs.blocks[i] for i = 1:Nblocks(rhs)])
end

Base.svd(A::AbstractBlock) = svd(A.block,thin=true)

# Default implementation

immutable VectorBlock <: AbstractVectorBlock
    block::Vector{Complex128}
end

immutable MatrixBlock <: AbstractMatrixBlock
    block::Matrix{Complex128}
end

immutable BlockVector <: AbstractBlockVector
    blocks::Vector{VectorBlock}
end

immutable BlockDiagonalMatrix <: AbstractBlockDiagonalMatrix
    blocks::Vector{MatrixBlock}
end

