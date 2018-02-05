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

struct AngularBlockVector <: BlockVector
    lmax :: Int
    m    :: Int
    blocks :: Vector{Vector{Complex128}}
end

function AngularBlockVector(lmax, m)
    blocks = Array{Vector{Complex128}}(lmax-m+1)
    AngularBlockVector(lmax, m, blocks)
end

Base.indices(vector::AngularBlockVector) = (vector.m:vector.lmax,)
Base.getindex(vector::AngularBlockVector, l) = vector.blocks[l-vector.m+1]
Base.setindex!(vector::AngularBlockVector, block, l) = vector.blocks[l-vector.m+1] = block

function Base.dot(lhs::AngularBlockVector, rhs::AngularBlockVector)
    sum(dot(b1, b2) for (b1, b2) in zip(lhs.blocks, rhs.blocks))
end

