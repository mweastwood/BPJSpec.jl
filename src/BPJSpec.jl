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

export GriddedVisibilities, grid!
export MModes, TransferMatrix
export tikhonov

using CasaCore.Measures
using CasaCore.Tables
using LibHealpix
using ProgressMeter
using TTCal

importall Base.Operators
import Cosmology
import GSL
import LibHealpix: Alm, lmax, mmax
import TTCal: Nfreq

include("special.jl")  # special functions
include("physics.jl")  # physical constants and cosmology
include("parallel.jl") # tools for parallel processing
#include("blocks.jl")  # block vectors and matrices
#include("itrf.jl")

# This function is useful to handle some of the special casing required for m == 0
two(m) = ifelse(m != 0, 2, 1)

# These functions define the filenames used when blocks are written to disk
block_filename(ν) = @sprintf("%.6fMHz.block", ν/1e6)
block_filename(m, ν) = @sprintf("%.6fMHz-%04d.block", ν/1e6, m)

# This function defines the ordering for blocks of a given m and frequency channel
# This ordering is assumed at various points so a careful check is needed before changing
block_index(mmax, m, channel) = (mmax+1)*(channel-1) + m + 1

include("visibilities.jl")
include("mmodes.jl")
include("transfermatrix.jl")
include("alm.jl")

#include("noise.jl")
#include("sky.jl")

end

