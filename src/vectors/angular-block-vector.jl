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
    mmax :: Int
    blocks :: Matrix{Vector{Complex128}}
end

function AngularBlockVector(lmax, mmax)
    blocks = Array{Vector{Complex128}}(lmax+1, mmax+1)
    AngularBlockVector(lmax, mmax, blocks)
end

function AngularBlockVector(input::SpectralBlockVector)
    lmax = mmax = input.mmax
    Nfreq = length(input.frequencies)
    output = AngularBlockVector(lmax, mmax)
    for m = 0:mmax, l = m:lmax
        output_block = zeros(Complex128, Nfreq)
        for β = 1:Nfreq
            input_block = input[m, β]
            output_block[β] = input_block[l-m+1]
        end
        output[l, m] = output_block
    end
    output
end

function SpectralBlockVector(input::AngularBlockVector)

end

indices(matrix::AngularBlockVector) =
    [(l, m) for m = 0:matrix.mmax for l = m:matrix.lmax]

Base.getindex(vector::AngularBlockVector, l, m) = vector.blocks[l+1, m+1]
Base.setindex!(vector::AngularBlockVector, block, l, m) = vector.blocks[l+1, m+1] = block

function Base.dot(lhs::AngularBlockVector, rhs::AngularBlockVector)
    output = complex(0.0)
    for m = 0:lhs.mmax, l = m:lhs.lmax
        output += dot(lhs[l, m], rhs[l, m])
    end
    output
end

