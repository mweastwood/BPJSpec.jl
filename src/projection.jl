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

immutable ProjectionMatrix{mmax}
    blocks::Vector{Matrix{Complex128}}
end

blocks(P::ProjectionMatrix) = P.blocks
block(P::ProjectionMatrix,m) = blocks(P)[abs(m)+1]
mmax{mmax}(P::ProjectionMatrix{mmax}) = mmax

@doc """
Construct a projection that projects the m-modes onto a
lower-dimensional space that retains as much information
about the sky as possible. This acts as a compression.

The heuristic I'm going to use here is that each block of
the matrix should be square. This will therefore be a
lossless compression in the sense that we're not losing
any information on the spherical harmonic coefficients,
because we're not actually discarding any of the singular
values.
""" ->
function ProjectionMatrix(transfermatrix::TransferMatrix)
    blocks = Array{Matrix{Complex128}}(mmax(transfermatrix)+1)
    for m = 0:mmax(transfermatrix)
        U,σ,V = svd(block(transfermatrix,m))
        blocks[m+1] = U[:,1:length(σ)]'
    end
    ProjectionMatrix{mmax(transfermatrix)}(blocks)
end

function *(P::ProjectionMatrix,B::TransferMatrix)
    if mmax(P) != mmax(B)
        error("The ProjectionMatrix and TransferMatrix must have the same mmax.")
    end
    Nbase′ = size(block(P,0),1)
    blocks = Array{Matrix{Complex128}}(mmax(B)+1)
    for m = 0:mmax(B)
        blocks[m+1] = block(P,m)*block(B,m)
    end
    TransferMatrix{Nbase′,lmax(B),mmax(B)}(blocks)
end

function *(P::ProjectionMatrix,v::MModes)
    if mmax(P) != mmax(v)
        error("The ProjectionMatrix and MModes must have the same mmax.")
    end
    Nbase′ = size(block(P,0),1)
    blocks = Array{Vector{Complex128}}(mmax(v)+1)
    for m = 0:mmax(v)
        blocks[m+1] = block(P,m)*block(v,m)
    end
    MModes{Nbase′,mmax(v)}(blocks)
end

