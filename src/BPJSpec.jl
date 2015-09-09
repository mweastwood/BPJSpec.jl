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

export ObsParam
export TransferMatrix
export MModes, visibilities, tikhonov
export ProjectionMatrix, compression
export CovarianceMatrix, ForegroundModel

importall Base.Operators
using HEALPix
import HEALPix: Alm, lmax, mmax
using CasaCore.Quanta
using CasaCore.Measures
using CasaCore.Tables
using MeasurementSets

include("special.jl")   # special functions
include("planewave.jl") # plane-wave expansion
include("physics.jl")   # physical constants and cosmology

# This function is useful to hand some of the
# special casing required for m == 0
two(m) = ifelse(m != 0, 2, 1)

include("obs.jl")
include("transfermatrix.jl")
include("mmodes.jl")
include("projection.jl")
include("alm.jl")
include("covariancematrix.jl")

include("cornerturn.jl")

end

