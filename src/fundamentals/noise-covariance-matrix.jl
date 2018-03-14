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

struct NoiseCovarianceMatrix{S} <: AbstractBlockMatrix{Diagonal{Float64}, 2}
    storage :: S
    cache   :: Cache{Diagonal{Float64}}
    mmax    :: Int
    frequencies :: Vector{typeof(1.0u"Hz")}
    bandwidth   :: Vector{typeof(1.0u"Hz")}
end
function NoiseCovarianceMatrix(storage::S, cache, mmax, frequencies, bandwidth) where S
    NoiseCovarianceMatrix{S}(storage, cache, mmax, frequencies, bandwidth)
end
metadata_fields(matrix::NoiseCovarianceMatrix) =
    (matrix.mmax, matrix.frequencies, matrix.bandwidth)
nblocks(::Type{<:NoiseCovarianceMatrix}, mmax, frequencies, bandwidth) =
    (mmax+1)*length(frequencies)
linear_index(matrix::NoiseCovarianceMatrix, m, β) =
    (matrix.mmax+1)*(β-1) + (m+1)
indices(matrix::NoiseCovarianceMatrix) =
    ((m, β) for β = 1:length(matrix.frequencies) for m = 0:matrix.mmax)

function compute!(matrix::NoiseCovarianceMatrix, model::NoiseModel, mmodes; progress=false)
    queue  = collect(indices(matrix))
    pool   = CachingPool(workers())
    if progress
        lck = ReentrantLock()
        prg = Progress(length(queue))
        increment() = (lock(lck); next!(prg); unlock(lck))
    end
    @sync for worker in workers()
        @async while length(queue) > 0
            m, β = shift!(queue)
            remotecall_fetch(_compute_noise_covariance, pool,
                             matrix, model, mmodes, m, β)
            progress && increment()
        end
    end
    matrix
end

function _compute_noise_covariance(matrix, model, mmodes, m, β)
    ν  = matrix.frequencies[β]
    Δν = matrix.bandwidth[β]
    σ  = ustrip(uconvert(u"Jy", model(m, ν, Δν)))
    N  = σ^2 .* ones(length(mmodes[m, β]))
    matrix[m, β] = Diagonal(N)
end

