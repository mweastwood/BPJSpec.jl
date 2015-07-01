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

module BPJSpec

export TransferMatrix, MModes
export blocks, block
export Nbase, Nfreq, lmax, mmax
export visibilities

using HEALPix
using CasaCore.Quanta
using CasaCore.Measures
using CasaCore.Tables

import HEALPix: Alm, lmax, mmax

const c = 2.99792e+8

include("special.jl") # special functions
include("planewave.jl") # plane-wave expansion

# These two functions are useful to
# handle some of the special casing
# required for m == 0 and negative m.
one(m) = ifelse(m  > 0, 1, 0)
two(m) = ifelse(m != 0, 2, 1)

include("transfermatrix.jl")
include("mmodes.jl")
include("alm.jl")

end

