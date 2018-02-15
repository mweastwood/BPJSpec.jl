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

#__precompile__()

module BPJSpec

# Matrices
export BlockDiagonalMatrix
export SpectralBlockDiagonalMatrix
export AngularCovarianceMatrix
export NoiseCovarianceMatrix
export HierarchicalTransferMatrix

# Vectors
export MModes

using CasaCore.Measures
using Cosmology
using Cubature
using FastTransforms
using FileIO, JLD2
using ProgressMeter
using StaticArrays
using Unitful, UnitfulAstro

T(A) = ctranspose(A)
H(A) = 0.5*(A+A') # guarantee Hermitian

"Try to make sure a matrix that should be positive definite is in fact positive definite."
function fix(A)
    N = size(A, 1)
    N == 0 && return A
    B = H(A)
    λ = eigvals(B)
    λmin = minimum(λ)
    λmax = maximum(λ)
    if λmin ≤ 0
        factor = N * eps(Float64) * λmax
        return B + factor * I
    else
        return B
    end
end

"Useful little function that helps account for grouping of positive and negative m."
two(m) = ifelse(m > 0, 2, 1)

include("parallel.jl")
include("cosmology.jl")
include("spherical-harmonics.jl")
include("metadata.jl")
include("hierarchy.jl")

abstract type SkyComponent end
include("sky/foregrounds.jl")
include("sky/signal.jl")
include("sky/noise.jl")

abstract type BlockMatrix end
include("matrices/block-diagonal-matrix.jl")
include("matrices/spectral-block-diagonal-matrix.jl")
include("matrices/angular-covariance-matrix.jl")
include("matrices/noise-covariance-matrix.jl")
include("matrices/transfer-matrix.jl")

abstract type BlockVector end
include("vectors/block-diagonal-vector.jl")
include("vectors/spectral-block-vector.jl")
include("vectors/angular-block-vector.jl")
include("vectors/random-angular-block-vector.jl")
include("vectors/white-noise-vector.jl")
include("vectors/random-vector.jl")

include("broadcasting.jl")

include("m-modes.jl")
#include("filter.jl")
#include("imaging.jl")
#include("fisher.jl")

include("algorithms/average-frequency-channels.jl")
include("algorithms/full-rank-compress.jl")

end

