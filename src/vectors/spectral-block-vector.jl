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

struct SpectralBlockVector <: BlockVector
    mmax :: Int
    frequencies :: Vector{typeof(1.0*u"Hz")}
    blocks :: Matrix{Vector{Complex128}}
end

function SpectralBlockVector(mmax, frequencies)
    blocks = Array{Vector{Complex128}}(mmax+1, length(frequencies))
    SpectralBlockVector(mmax, frequencies, blocks)
end

Base.indices(vector::SpectralBlockVector) = (0:vector.mmax, 1:length(vector.frequencies))
Base.getindex(vector::SpectralBlockVector, m, β) = vector.blocks[m+1, β]
Base.setindex!(vector::SpectralBlockVector, block, m, β) = vector.blocks[m+1, β] = block

function Base.getindex(vector::SpectralBlockVector, m)
    blocks = [vector[m, ν] for ν in vector.frequencies]
    X = sum(length.(blocks))
    output = zeros(Complex128, X)
    x = 1
    for block in blocks
        output[x:x+size(block, 1)-1] = block
        x += size(block, 1)
    end
    output
end

