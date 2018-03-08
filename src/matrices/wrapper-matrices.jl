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

Base.show(io::IO, matrix::AbstractBlockMatrix) = show(io, matrix.matrix)
indices(matrix::AbstractBlockMatrix) = indices(matrix.matrix.metadata)

cache!(matrix::AbstractBlockMatrix) = cache!(matrix.matrix)
flush!(matrix::AbstractBlockMatrix) = flush!(matrix.matrix)

distribute_write(matrix::AbstractBlockMatrix) = distribute_write(matrix.matrix)
distribute_read(matrix::AbstractBlockMatrix) = distribute_read(matrix.matrix)

storage(matrix::AbstractBlockMatrix) = matrix.matrix.storage
frequencies(matrix::AbstractBlockMatrix) = matrix.matrix.metadata.frequencies
bandwidth(matrix::AbstractBlockMatrix) = matrix.matrix.metadata.bandwidth
position(matrix::AbstractBlockMatrix) = matrix.matrix.metadata.position
baselines(matrix::AbstractBlockMatrix) = matrix.matrix.metadata.baselines
phase_center(matrix::AbstractBlockMatrix) = matrix.matrix.metadata.phase_center
lmax(matrix::AbstractBlockMatrix) = matrix.matrix.metadata.lmax
mmax(matrix::AbstractBlockMatrix) = matrix.matrix.metadata.mmax

"Metadata for matrices that will be split into blocks of m."
struct MMax <: MatrixMetadata
    mmax :: Int
end
Base.show(io::IO, metadata::MMax) = @printf(io, "{mmax: %d}", metadata.mmax)
indices(metadata::MMax) = 0:metadata.mmax
number_of_blocks(metadata::MMax) = metadata.mmax+1

"Metadata for matrices that will be split into blocks of m and frequency."
struct MMaxFrequencies <: MatrixMetadata
    mmax        :: Int
    frequencies :: Vector{typeof(1.0*u"Hz")}
    bandwidth   :: Vector{typeof(1.0*u"Hz")}
end
function Base.show(io::IO, metadata::MMaxFrequencies)
    @printf(io, "{mmax: %d, ν: %.3f MHz..%.3f MHz, Δν: %.3f MHz total}",
            metadata.mmax,
            ustrip(uconvert(u"MHz", metadata.frequencies[1])),
            ustrip(uconvert(u"MHz", metadata.frequencies[end])),
            ustrip(uconvert(u"MHz", sum(metadata.bandwidth))))
end
function indices(metadata::MMaxFrequencies)
    ((m, β) for β = 1:length(metadata.frequencies) for m = 0:metadata.mmax)
end
function number_of_blocks(metadata::MMaxFrequencies)
    (metadata.mmax+1, length(metadata.frequencies))
end

"Matrix that is split into blocks of m."
struct MBlockMatrix{S} <: AbstractBlockMatrix
    matrix :: BlockMatrix{Matrix{Complex128}, 1, MMax, S}
end
function MBlockMatrix(storage::Mechanism, mmax)
    metadata = MMax(mmax)
    matrix = BlockMatrix{Matrix{Complex128}, 1}(storage, metadata)
    MBlockMatrix(matrix)
end
MBlockMatrix(path::String) = MBlockMatrix(BlockMatrix{Matrix{Complex128}, 1}(path))

Base.getindex(matrix::MBlockMatrix, m) = matrix.matrix[m+1]
Base.setindex!(matrix::MBlockMatrix, block, m) = matrix.matrix[m+1] = block

"Matrix that is split into blocks of m and frequency."
struct MFBlockMatrix{S} <: AbstractBlockMatrix
    matrix :: BlockMatrix{Matrix{Complex128}, 2, MMaxFrequencies, S}
end
function MFBlockMatrix(storage::Mechanism, mmax, frequencies, bandwidth)
    metadata = MMaxFrequencies(mmax, frequencies, bandwidth)
    matrix = BlockMatrix{Matrix{Complex128}, 2}(storage, metadata)
    MFBlockMatrix(matrix)
end
MFBlockMatrix(path::String) = MFBlockMatrix(BlockMatrix{Matrix{Complex128}, 2}(path))

Base.getindex(matrix::MFBlockMatrix, m, β) = matrix.matrix[m+1, β]
Base.setindex!(matrix::MFBlockMatrix, block, m, β) = matrix.matrix[m+1, β] = block

"Vector that is split into blocks of m."
struct MBlockVector{S} <: AbstractBlockMatrix
    matrix :: BlockMatrix{Vector{Complex128}, 1, MMax, S}
end
function MBlockVector(storage::Mechanism, mmax)
    metadata = MMax(mmax)
    matrix = BlockMatrix{Vector{Complex128}, 1}(storage, metadata)
    MBlockVector(matrix)
end
MBlockVector(path::String) = MBlockVector(BlockMatrix{Vector{Complex128}, 1}(path))

Base.getindex(vector::MBlockVector, m) = vector.matrix[m+1]
Base.setindex!(vector::MBlockVector, block, m) = vector.matrix[m+1] = block

"Vector that is split into blocks of m and frequency."
struct MFBlockVector{S} <: AbstractBlockMatrix
    matrix :: BlockMatrix{Vector{Complex128}, 2, MMaxFrequencies, S}
end
function MFBlockVector(storage::Mechanism, mmax, frequencies, bandwidth)
    metadata = MMaxFrequencies(mmax, frequencies, bandwidth)
    matrix = BlockMatrix{Vector{Complex128}, 2}(storage, metadata)
    MFBlockVector(matrix)
end
MFBlockVector(path::String) = MFBlockVector(BlockMatrix{Vector{Complex128}, 2}(path))

Base.getindex(vector::MFBlockVector, m, β) = vector.matrix[m+1, β]
Base.setindex!(vector::MFBlockVector, block, m, β) = vector.matrix[m+1, β] = block

