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

const NoiseCovarianceMatrix = SpectralBlockDiagonalMatrix{Diagonal{Float64}}

function NoiseCovarianceMatrix(path, mmax, frequencies, bandwidth, hierarchy, model::NoiseModel)
    output = NoiseCovarianceMatrix(path, mmax, frequencies, bandwidth)
    compute!(output, hierarchy, model)
    output
end

function compute!(matrix::NoiseCovarianceMatrix, hierarchy, model::NoiseModel)
    Nfreq = length(matrix.frequencies)
    for β = 1:Nfreq
        ν  = matrix.frequencies[β]
        Δν = matrix.bandwidth[β]
        for m = 0:matrix.mmax
            σ = ustrip(uconvert(u"Jy", model(m, ν, Δν)))
            N = σ^2 .* ones(two(m)*Nbase(hierarchy, m))
            matrix[m, β] = Diagonal(N)
        end
    end
end

Base.show(io::IO, matrix::NoiseCovarianceMatrix) =
    print(io, "NoiseCovarianceMatrix: ", matrix.path)

