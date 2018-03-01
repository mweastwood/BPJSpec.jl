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

doc"""
    module BPJSpec

BPJSpec is a 21-cm power spectrum code developed for the OVRO-LWA based on the $m$-mode analysis
formalism.

## $m$-Mode Analysis

$m$-mode analysis is a relatively new technique that has been developed for drift-scanning
telescopes. This technique begins by applying a Fourier transform over sidereal time to the measured
visibilities. This transformation introduces block-diagonal structure into the transfer matrix,
which describes how the interferometer responds to the sky. The additional sparsity in the transfer
matrix facilitates a matrix-algebra approach to interferometric imaging and power spectrum
estimation.

## References

* https://arxiv.org/abs/1302.0327
* https://arxiv.org/abs/1401.2095
* https://arxiv.org/abs/1711.00466
"""
module BPJSpec

# Matrices
export BlockDiagonalMatrix, DenseBlockDiagonalMatrix
export SpectralBlockDiagonalMatrix, DenseSpectralBlockDiagonalMatrix
export AngularCovarianceMatrix, NoiseCovarianceMatrix
export TransferMatrix, HierarchicalTransferMatrix

# Vectors
export SpectralBlockVector
export MModes

using Unitful, UnitfulAstro # Travis CI fails with "invalid age range update" unless this is first

using ApproxFun
using CasaCore.Measures
using Cosmology
using Cubature
using FastTransforms
using FileIO, JLD2
using ProgressMeter
using StaticArrays

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
two(m) = ifelse(m != 0, 2, 1)

include("parallel.jl")
include("cosmology.jl")
include("spherical-harmonics.jl")
include("metadata.jl")
include("hierarchy.jl")

abstract type SkyComponent end
struct NoComponent <: SkyComponent end
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
include("algorithms/average-frequency-channels.jl")
include("algorithms/full-rank-compress.jl")
include("algorithms/karhunen-loeve-transforms.jl")

include("m-modes.jl")
include("imaging.jl")
include("fisher.jl")

end

