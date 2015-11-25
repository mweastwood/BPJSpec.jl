# Copyright (c) 2015 Michael Eastwood
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

__precompile__()

module BPJSpec

export itrf_baselines, itrf_phasecenter, itrf_beam
export create_empty_visibilities, grid_visibilities, load_visibilities

export visibilities, mmodes, transfer, tikhonov
export preserve_singular_values

export NoiseModel, ForegroundModel, SignalModel
export covariance_matrix, dimensionful_powerspectrum

using CasaCore.Measures
using CasaCore.Tables
using HDF5, JLD
using LibHealpix
using ProgressMeter
using TTCal

importall Base.Operators
import Cosmology
import GSL
import LibHealpix: Alm, lmax, mmax

include("special.jl") # special functions
include("physics.jl") # physical constants and cosmology
include("blocks.jl")  # block vectors and matrices
include("itrf.jl")

# This function is useful to handle some of the
# special casing required for m == 0
two(m) = ifelse(m != 0, 2, 1)

include("visibilities.jl")
include("transfermatrix.jl")
include("mmodes.jl")
include("alm.jl")
include("noise.jl")

include("covariancematrix.jl")

function load(filename, args...)
    local description
    jldopen(filename,"r") do file
        description = read(file["description"])
    end
    if description == "transfer matrix"
        return load(filename, TransferMeta(args...))
    elseif description == "m-modes"
        return load(filename, MModesMeta(args...))
    elseif description == "noise covariance matrix"
        return load(filename, NoiseMeta(args...))
    else
        error("Unrecognized format.")
    end
end

end

