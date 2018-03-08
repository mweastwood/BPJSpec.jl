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

struct NoiseCovarianceMatrix{S} <: AbstractBlockMatrix
    matrix :: BlockMatrix{Diagonal{Float64}, 2, MMaxFrequencies, S}
end
function NoiseCovarianceMatrix(storage::Mechanism, mmax, frequencies, bandwidth)
    metadata = MMaxFrequencies(mmax, frequencies, bandwidth)
    matrix = BlockMatrix{Diagonal{Float64}, 2}(storage, metadata)
    NoiseCovarianceMatrix(matrix)
end
function NoiseCovarianceMatrix(path::String)
    NoiseCovarianceMatrix(BlockMatrix{Diagonal{Float64}, 2}(path))
end

Base.getindex(matrix::NoiseCovarianceMatrix, m, β) = matrix.matrix[m+1, β]
Base.setindex!(matrix::NoiseCovarianceMatrix, block, m, β) = matrix.matrix[m+1, β] = block

function compute!(matrix::NoiseCovarianceMatrix, hierarchy::Hierarchy, model::NoiseModel)
    Nfreq = length(frequencies(matrix))
    for β = 1:Nfreq
        ν  = frequencies(matrix)[β]
        Δν = bandwidth(matrix)[β]
        for m = 0:mmax(matrix)
            σ = ustrip(uconvert(u"Jy", model(m, ν, Δν)))
            N = σ^2 .* ones(two(m)*Nbase(hierarchy, m))
            matrix[m, β] = Diagonal(N)
        end
    end
end

