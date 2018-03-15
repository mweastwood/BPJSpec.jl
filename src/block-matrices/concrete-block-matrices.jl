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

"Array that is split into arbitrary blocks."
struct SimpleBlockArray{T, N, S} <: AbstractBlockMatrix{Array{T, N}, 1}
    storage :: S
    cache   :: Cache{Array{T, N}}
    length  :: Int
end
metadata_fields(array::SimpleBlockArray) = (array.length,)
nblocks(::Type{<:SimpleBlockArray}, length) = length
linear_index(::SimpleBlockArray, idx) = idx
indices(array::SimpleBlockArray) = 1:array.length

"Array that is split into blocks of m."
struct MBlockArray{T, N, S} <: AbstractBlockMatrix{Array{T, N}, 1}
    storage :: S
    cache   :: Cache{Array{T, N}}
    mmax    :: Int
end
metadata_fields(array::MBlockArray) = (array.mmax,)
nblocks(::Type{<:MBlockArray}, mmax) = mmax+1
linear_index(::MBlockArray, m) = m+1
indices(array::MBlockArray) = 0:array.mmax

"Array that is split into blocks of frequency."
struct FBlockArray{T, N, S} <: AbstractBlockMatrix{Array{T, N}, 1}
    storage :: S
    cache   :: Cache{Array{T, N}}
    frequencies :: Vector{typeof(1.0u"Hz")}
    bandwidth   :: Vector{typeof(1.0u"Hz")}
end
metadata_fields(array::FBlockArray) = (array.frequencies, array.bandwidth)
nblocks(::Type{<:FBlockArray}, frequencies, bandwidth) = length(frequencies)
linear_index(::FBlockArray, β) = β
indices(array::FBlockArray) = 1:length(array.frequencies)

"Array that is split into blocks of m and frequency."
struct MFBlockArray{T, N, S} <: AbstractBlockMatrix{Array{T, N}, 2}
    storage :: S
    cache   :: Cache{Array{T, N}}
    mmax    :: Int
    frequencies :: Vector{typeof(1.0u"Hz")}
    bandwidth   :: Vector{typeof(1.0u"Hz")}
end
metadata_fields(array::MFBlockArray) = (array.mmax, array.frequencies, array.bandwidth)
nblocks(::Type{<:MFBlockArray}, mmax, frequencies, bandwidth) = (mmax+1)*length(frequencies)
linear_index(array::MFBlockArray, m, β) = (array.mmax+1)*(β-1) + (m+1)
indices(array::MFBlockArray) = ((m, β) for β = 1:length(array.frequencies) for m = 0:array.mmax)

"Array that is split into blocks of l."
struct LBlockArray{T, N, S} <: AbstractBlockMatrix{Array{T, N}, 1}
    storage :: S
    cache   :: Cache{Array{T, N}}
    lmax    :: Int
    frequencies :: Vector{typeof(1.0u"Hz")}
    bandwidth   :: Vector{typeof(1.0u"Hz")}
end
metadata_fields(array::LBlockArray) = (array.lmax, array.frequencies, array.bandwidth)
nblocks(::Type{<:LBlockArray}, lmax, frequencies, bandwidth) = lmax+1
linear_index(array::LBlockArray, l) = l+1
indices(array::LBlockArray) = L(0):L(array.lmax)

"Array that is split into blocks of l and m."
struct LMBlockArray{T, N, S} <: AbstractBlockMatrix{Array{T, N}, 2}
    storage :: S
    cache   :: Cache{Array{T, N}}
    lmax    :: Int
    mmax    :: Int
    frequencies :: Vector{typeof(1.0u"Hz")}
    bandwidth   :: Vector{typeof(1.0u"Hz")}
end
metadata_fields(array::LMBlockArray) =
    (array.lmax, array.mmax, array.frequencies, array.bandwidth)
nblocks(::Type{<:LMBlockArray}, lmax, mmax, frequencies, bandwidth) =
    ((2lmax + 2 - mmax) * (mmax + 1)) ÷ 2
linear_index(array::LMBlockArray, l, m) =
    (m * (2array.lmax - m + 3)) ÷ 2 + l - m + 1
indices(array::LMBlockArray) =
    ((l, m) for m = 0:array.mmax for l = m:array.lmax)

const SimpleBlockVector = SimpleBlockArray{Complex128, 1}
const SimpleBlockMatrix = SimpleBlockArray{Complex128, 2}
const  MBlockVector =  MBlockArray{Complex128, 1}
const  MBlockMatrix =  MBlockArray{Complex128, 2}
const  FBlockVector =  FBlockArray{Complex128, 1}
const  FBlockMatrix =  FBlockArray{Complex128, 2}
const MFBlockVector = MFBlockArray{Complex128, 1}
const MFBlockMatrix = MFBlockArray{Complex128, 2}
# The following type uses Matrix{Float64} blocks because we want to use it as an angular covariance
# matrix, which is block diagonal in l and has real elements.
const  LBlockMatrix =  LBlockArray{   Float64, 2}
const LMBlockVector = LMBlockArray{Complex128, 1}

# Specialized indexing rules
# ==========================

function Base.getindex(matrix::MFBlockVector, m::Int)
    blocks = [matrix[m, β] for β = 1:length(matrix.frequencies)]
    X = sum(size.(blocks, 1))
    output = zeros(eltype(first(blocks)), X)
    x = 1
    for block in blocks
        output[x:x+size(block, 1)-1] = block
        x += size(block, 1)
    end
    output
end

function Base.getindex(matrix::MFBlockMatrix, m::Int)
    blocks = [matrix[m, β] for β = 1:length(matrix.frequencies)]
    X = sum(size.(blocks, 1))
    Y = sum(size.(blocks, 2))
    output = zeros(eltype(first(blocks)), X, Y)
    x = y = 1
    for block in blocks
        output[x:x+size(block, 1)-1, y:y+size(block, 2)-1] = block
        x += size(block, 1)
        y += size(block, 2)
    end
    output
end

