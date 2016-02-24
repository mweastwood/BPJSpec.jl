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

abstract Block
abstract Metadata

immutable NoMetadata <: Metadata end

"""
    immutable VectorBlock <: Block

A block that composes a vector.
"""
immutable VectorBlock <: Block
    block::Vector{Complex128}
end

VectorBlock(sz) = VectorBlock(zeros(Complex128,sz))

"""
    immutable MatrixBlock <: Block

A dense block that composes a matrix.
"""
immutable MatrixBlock <: Block
    block::Matrix{Complex128}
end

MatrixBlock(sz) = MatrixBlock(zeros(Complex128,sz))

"""
    immutable DiagonalMatrixBlock <: Block

A diagonal block that composes a matrix.
"""
immutable DiagonalMatrixBlock <: Block
    block::Diagonal{Complex128}
end

DiagonalMatrixBlock(sz) = DiagonalMatrixBlock(Diagonal(zeros(Complex128,sz)))

"""
    immutable Blocks{T<:Block, S<:Metadata}

A list of blocks. This can represent a block diagonal matrix,
or a vector that is organized into blocks.
"""
immutable Blocks{T<:Block, S<:Metadata}
    blocks::Vector{T}
    meta::S
end

Blocks{T}(blocks::Vector{T}) = Blocks{T,NoMetadata}(blocks,NoMetadata())

typealias BlockVector{S} Blocks{VectorBlock, S}
typealias BlockMatrix{S} Union{Blocks{MatrixBlock, S}, Blocks{DiagonalMatrixBlock, S}}

getindex(A::VectorBlock,i) = A.block[i]
setindex!(A::VectorBlock,x,i) = A.block[i] = x

getindex(A::MatrixBlock,i,j) = A.block[i,j]
setindex!(A::MatrixBlock,x,i,j) = A.block[i,j] = x

getindex(A::DiagonalMatrixBlock,i,j) = A.block[i,j]
setindex!(A::DiagonalMatrixBlock,x,i,j) = A.block[i,j] = x

getindex(A::Blocks,i) = A.blocks[i]

==(A::Block,B::Block) = A.block == B.block
==(A::Blocks,B::Blocks) = A.blocks == B.blocks && A.meta == B.meta

Base.length(A::Block) = length(A.block)
Base.size(A::Block) = size(A.block)

function Base.length(A::BlockVector)
    sz = 0
    for block in A.blocks
        sz += length(block)
    end
    sz
end

function Base.size(A::BlockMatrix)
    sz1 = sz2 = 0
    for block in A.blocks
        sz = size(block)
        sz1 += sz[1]
        sz2 += sz[2]
    end
    sz1,sz2
end

function Base.full(A::BlockVector)
    B = zeros(Complex128,length(A))
    idx = 1
    for block in A.blocks
        sz = length(block)
        B[idx:idx+sz-1] = block.block
        idx += sz
    end
    B
end

function Base.full(A::BlockMatrix)
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

ctranspose{T<:Block}(A::T) = T(A.block')
ctranspose(A::Blocks) = Blocks([A.blocks[i]' for i = 1:length(A.blocks)],A.meta)

*(A::MatrixBlock, B::DiagonalMatrixBlock) = MatrixBlock(A.block * B.block)
*(A::DiagonalMatrixBlock, B::MatrixBlock) = MatrixBlock(A.block * B.block)
*(A::MatrixBlock, B::DiagonalMatrixBlock, C::MatrixBlock) = MatrixBlock(A.block * B.block * C.block)
*(A::MatrixBlock, B::MatrixBlock) = MatrixBlock(A.block * B.block)
*(A::MatrixBlock, B::VectorBlock) = VectorBlock(A.block * B.block)

function *(A::Blocks, B::Blocks)
    Blocks([A.blocks[i]*B.blocks[i] for i = 1:length(B.blocks)],B.meta)
end

function *(A::Blocks, B::Blocks, C::Blocks)
    Blocks([A.blocks[i]*B.blocks[i]*C.blocks[i] for i = 1:length(B.blocks)],B.meta)
end

Base.svd(A::Block) = svd(A.block,thin=true)

