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

struct NoiseCovarianceMatrix end

function create(::Type{NoiseCovarianceMatrix}, path::String,
                model::NoiseModel, metadata::Metadata, hierarchy::Hierarchy;
                rm=false, progress=false)
    storage = MultipleFiles(path)
    output  = create(MFDiagonalBlockMatrix, storage, hierarchy.lmax,
                     metadata.frequencies, metadata.bandwidth, rm=rm)
    compute!(NoiseCovarianceMatrix, output, model, hierarchy, progress=progress)
    output
end

function compute!(::Type{NoiseCovarianceMatrix}, matrix::MFDiagonalBlockMatrix,
                  model::NoiseModel, hierarchy::Hierarchy; progress=false)
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
                             matrix, model, hierarchy, m, β)
            progress && increment()
        end
    end
    matrix
end

function _compute_noise_covariance(matrix, model, hierarchy, m, β)
    ν  = matrix.frequencies[β]
    Δν = matrix.bandwidth[β]
    σ  = ustrip(uconvert(u"Jy", model(m, ν, Δν)))
    N  = σ^2 .* ones(two(m)*Nbase(hierarchy, m))
    matrix[m, β] = Diagonal(N)
end

