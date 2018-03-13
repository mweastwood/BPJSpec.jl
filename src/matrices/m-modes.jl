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

doc"""
    struct MModes{S} <: AbstractBlockMatrix{Vector{Complex128}, 2}

This type represents the $m$-modes measured by the interferometer.

# Fields

* `storage` contains instructions on how to read the m-modes from disk
* `cache` is used if we want to keep the $m$-modes in memory
* `mmax` is the largest value of the $m$ quantum number
* `frequencies` is the list of frequencies
* `bandwidth` is the bandwidth associated with each frequency channel
"""
struct MModes{S} <: AbstractBlockMatrix{Vector{Complex128}, 2}
    storage :: S
    cache   :: Cache{Vector{Complex128}}
    mmax    :: Int
    frequencies :: Vector{typeof(1.0u"Hz")}
    bandwidth   :: Vector{typeof(1.0u"Hz")}
end
function MModes(storage::S, cache, mmax, frequencies, bandwidth) where S
    MModes{S}(storage, cache, mmax, frequencies, bandwidth)
end
metadata_fields(mmodes::MModes) = (mmodes.mmax, mmodes.frequencies, mmodes.bandwidth)
nblocks(::Type{<:MModes}, mmax, frequencies, bandwidth) = (mmax+1)*length(frequencies)
linear_index(mmodes::MModes, m, β) = (mmodes.mmax+1)*(β-1) + (m+1)
indices(mmodes::MModes) = ((m, β) for β = 1:length(mmodes.frequencies) for m = 0:mmodes.mmax)

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
    for m = 1:mmodes.mmax
        for α = 1:Nbase
            α1 = 2α-1 # positive m
            α2 = 2α-0 # negative m
            block[α1] =      transformed_visibilities[      m+1, α]
            block[α2] = conj(transformed_visibilities[Ntime+1-m, α])
        end
        mmodes[m, β] = block
    end
end

