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

# m-Mode Analysis

m-mode analysis is a relatively new technique that has been developed for drift-scanning telescopes.
This technique begins by applying a Fourier transform over sidereal time to the measured
visibilities. This transformation introduces block-diagonal structure into the transfer matrix,
which describes how the interferometer responds to the sky. The additional sparsity in the transfer
matrix facilitates a matrix-algebra approach to interferometric imaging and power spectrum
estimation.

# References

* https://arxiv.org/abs/1302.0327
* https://arxiv.org/abs/1401.2095
* https://arxiv.org/abs/1711.00466
"""
module BPJSpec

export NoFile, SingleFile, MultipleFiles
export SimpleBlockMatrix, SimpleBlockVector
export MBlockMatrix, FBlockMatrix, MFBlockMatrix, LBlockMatrix
export MBlockVector, FBlockVector, MFBlockVector, LMBlockVector
export TransferMatrix, NoiseCovarianceMatrix, AngularCovarianceMatrix
export MModes, MultiFrequencyAlm
export RandomBlockVector, WhiteNoiseBlockVector
export ProgressBar, create, compute!, cache!, flush!

using Unitful, UnitfulAstro # Travis CI fails with "invalid age range update" unless this is first

using ApproxFun
using CasaCore.Measures
using Cosmology
using Cubature
using FastTransforms
using JLD2
using MacroTools
using ProgressMeter
using StaticArrays

# Defines an extension of FastTransforms.jl that provides a more convenient interface for fast
# spherical harmonic transforms.
include("wrappers/FastTransformsWrapper.jl")
using .FastTransformsWrapper

# Defines an interface to Cosmology.jl that attaches units from UnitfulAstro.jl. Furthermore, we
# include functionality from ApproxFun.jl to allow more rapid evaluation of these functions.
include("wrappers/CosmologyWrapper.jl")
using .CosmologyWrapper

# Defines an itnerface to GSL.jl's spherical harmonic functions. GSL is a pretty large dependency,
# but it provides robust spherical harmonic routines that I am confident I can rely on if necessary.
include("wrappers/GSLWrapper.jl")
using .GSLWrapper

include("utilities/misc.jl")
include("utilities/parallel.jl")
include("utilities/recombination-lines.jl")
include("utilities/white-noise.jl")

abstract type SkyComponent end
include("sky-components/foregrounds.jl")
include("sky-components/signal.jl")

include("interferometer/metadata.jl")
include("interferometer/baseline-hierarchy.jl")
include("interferometer/noise-model.jl")

include("block-matrices/storage-mechanisms.jl")
include("block-matrices/abstract-block-matrix.jl")
include("block-matrices/concrete-block-matrices.jl")
include("block-matrices/broadcasting.jl")

include("fundamentals/transfer-matrix.jl")
include("fundamentals/noise-covariance-matrix.jl")
include("fundamentals/angular-covariance-matrix.jl")
include("fundamentals/m-modes.jl")
include("fundamentals/multi-frequency-alm.jl")
include("fundamentals/random-block-vector.jl")

include("algorithms/permute-m-modes.jl")
include("algorithms/propagate-flags.jl")
include("algorithms/average-frequency-channels.jl")
include("algorithms/full-rank-compress.jl")
include("algorithms/karhunen-loeve-transforms.jl")
include("algorithms/tikhonov-regularization.jl")

include("quadratic-estimator/fisher-information.jl")
include("quadratic-estimator/noise-bias.jl")
include("quadratic-estimator/q-estimator.jl")
include("quadratic-estimator/mixing-matrix.jl")

end

