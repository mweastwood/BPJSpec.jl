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

"""
    ProjectionMatrixBlock{mmax}

A single block of a `ProjectionMatrix{mmax}` corresponding to
a single value of `m`.
"""
immutable ProjectionMatrixBlock{mmax}
    m::Int
    block::Matrix{Complex128}
end

"""
    ProjectionMatrix{mmax}

This type represents a block-diagonal matrix projection.
It is used to compress the dataset and remove unwanted
foreground power.
"""
immutable ProjectionMatrix{mmax}
    blocks::Vector{ProjectionMatrixBlock{mmax}}
end

mmax{m}(P::ProjectionMatrixBlock{m}) = m
mmax{m}(P::ProjectionMatrix{m}) = m
getindex(P::ProjectionMatrix,m) = P.blocks[m+1]

function *(P::ProjectionMatrixBlock,B::TransferMatrixBlock)
    mmax(P) == mmax(B) || error("The values of mmax must be the same.")
    P.m == B.m || error("The values of m must be the same.")
    TransferMatrixBlock{lmax(B),mmax(B)}(B.m,P.block*B.block)
end

function *(P::ProjectionMatrix,B::TransferMatrix)
    mmax(P) == mmax(B) || error("The values of mmax must be the same.")
    blocks = [P[m]*B[m] for m = 0:mmax(B)]
    TransferMatrix{lmax(B),mmax(B)}(blocks)
end

function *(P::ProjectionMatrixBlock,v::MModeBlock)
    mmax(P) == mmax(v) || error("The values of mmax must be the same.")
    P.m == v.m || error("The values of m must be the same.")
    MModeBlock{mmax(v)}(v.m,P.block*v.block)
end

function *(P::ProjectionMatrix,v::MModes)
    mmax(P) == mmax(v) || error("The values of mmax must be the same.")
    blocks = [P[m]*v[m] for m = 0:mmax(v)]
    MModes{mmax(v)}(blocks)
end

"""
    compression(B::TransferMatrix) -> ProjectionMatrix

Construct a projection that projects the m-modes onto a
lower-dimensional space that retains as much information
about the sky as possible. This acts as a compression.

The heuristic I'm going to use here is that each block of
the matrix should be square. This will therefore be a
lossless compression in the sense that we're not losing
any information on the spherical harmonic coefficients,
because we're not actually discarding any of the singular
values.
"""
function compression(B::TransferMatrix)
    blocks = Array{ProjectionMatrixBlock{mmax(B)}}(mmax(B)+1)
    for m = 0:mmax(B)
        U,σ,V = svd(B[m])
        blocks[m+1] = ProjectionMatrixBlock{mmax(B)}(m,U[:,1:length(σ)]')
    end
    ProjectionMatrix{mmax(B)}(blocks)
end

