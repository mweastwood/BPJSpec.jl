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

"""
    struct MModes{S} <: AbstractBlockMatrix

This type represents the m-modes measured by the interferometer.

# Fields

* `matrix` is the `BlockMatrix` backend that stores the m-modes
"""
struct MModes{S} <: AbstractBlockMatrix
    matrix :: BlockMatrix{Vector{Complex128}, 2, MMaxFrequencies, S}
end
function MModes(storage::Mechanism, mmax, frequencies, bandwidth)
    metadata = MMaxFrequencies(mmax, frequencies, bandwidth)
    matrix = BlockMatrix{Vector{Complex128}, 2}(storage, metadata)
    MModes(matrix)
end
MModes(path::String) = MModes(BlockMatrix{Vector{Complex128}, 2}(path))

Base.getindex(vector::MModes, m, β) = vector.matrix[m+1, β]
Base.setindex!(vector::MModes, block, m, β) = vector.matrix[m+1, β] = block

"Compute m-modes from two dimensional matrix of visibilities (time × baseline)."
function compute!(mmodes::MModes, visibilities::Matrix{Complex128}, β)
    store!(mmodes, fourier_transform(visibilities), β)
end

function fourier_transform(matrix)
    Ntime, Nbase = size(matrix)
    planned_fft = plan_fft(matrix, 1)
    (planned_fft*matrix) ./ Ntime
end

function store!(mmodes, transformed_visibilities, β)
    Ntime, Nbase = size(transformed_visibilities)

    # m = 0
    block = zeros(Complex128, Nbase)
    for α = 1:Nbase
        block[α] = transformed_visibilities[1, α]
    end
    mmodes[0, β] = block

    # m > 0
    block = zeros(Complex128, 2Nbase)
    for m = 1:mmax(mmodes)
        for α = 1:Nbase
            α1 = 2α-1 # positive m
            α2 = 2α-0 # negative m
            block[α1] =      transformed_visibilities[      m+1, α]
            block[α2] = conj(transformed_visibilities[Ntime+1-m, α])
        end
        mmodes[m, β] = block
    end
end

