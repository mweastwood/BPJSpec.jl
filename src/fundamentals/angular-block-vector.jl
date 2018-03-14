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

struct AngularBlockVector{S} <: AbstractBlockMatrix{Vector{Complex128}, 2}
    storage :: S
    cache   :: Cache{Vector{Complex128}}
    lmax :: Int
    mmax :: Int
end
function AngularBlockVector(storage::S, cache, lmax, mmax) where S
    AngularBlockVector{S}(storage, cache, lmax, mmax)
end
metadata_fields(matrix::AngularBlockVector) = (matrix.lmax, matrix.mmax)
nblocks(::Type{<:AngularBlockVector}, lmax, mmax) = ((2lmax + 2 - mmax) * (mmax + 1)) ÷ 2
linear_index(matrix::AngularBlockVector, l, m) = (m * (2matrix.lmax - m + 3)) ÷ 2 + l - m + 1
indices(matrix::AngularBlockVector) = ((l, m) for m = 0:matrix.mmax for l = m:matrix.lmax)

#function AngularBlockVector(input::SpectralBlockVector)
#    lmax = mmax = input.mmax
#    Nfreq = length(input.frequencies)
#    output = AngularBlockVector(lmax, mmax)
#    for m = 0:mmax, l = m:lmax
#        output_block = zeros(Complex128, Nfreq)
#        for β = 1:Nfreq
#            input_block = input[m, β]
#            output_block[β] = input_block[l-m+1]
#        end
#        output[l, m] = output_block
#    end
#    output
#end

function AngularBlockVector(input::MBlockVector)
    lmax = mmax = input.mmax
    Nfreq = length(input[0]) ÷ (lmax+1)
    output = create(AngularBlockVector, lmax, mmax)
    for m = 0:mmax, l = m:lmax
        output_block = zeros(Complex128, Nfreq)
        for β = 1:Nfreq
            input_block = input[m]
            output_block[β] = input_block[(lmax-m+1)*(β-1) + (l-m+1)]
        end
        output[l, m] = output_block
    end
    output
end

function Base.dot(lhs::AngularBlockVector, rhs::AngularBlockVector)
    output = zero(Complex128)
    for m = 0:lhs.mmax, l = m:lhs.lmax
        output += dot(lhs[l, m], rhs[l, m])
    end
    output
end

