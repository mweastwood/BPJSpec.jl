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

#include("spherical-harmonics.jl")

#using JLD2
#using Unitful, UnitfulAstro

using JLD2
using StaticArrays
using LibHealpix

struct SphericalHarmonicMetadata
    lmax :: Int
    mmax :: Int
    ν :: Vector{typeof(1.0*u"Hz")}
    function SphericalHarmonicMetadata(lmax, mmax, ν)
        mmax ≤ lmax || throw(ArgumentError("spherical harmonics require mmax ≤ lmax"))
        new(lmax, mmax, ν)
    end
end

struct InterferometerMetadata
    u :: Vector{typeof{1.0*u"m"}}
    v :: Vector{typeof{1.0*u"m"}}
    w :: Vector{typeof{1.0*u"m"}}
    beam :: RingHealpixMap{Float64}
    phase_center :: SVector{3, Float64}
end

include("transfermatrix.jl")

#export GriddedVisibilities, grid!
#export MModes, TransferMatrix
#export tikhonov
#
#using CasaCore.Measures
#using CasaCore.Tables
#using LibHealpix
#using ProgressMeter
#using TTCal
#
#importall Base.Operators
#import Cosmology
#import GSL
#import LibHealpix: Alm, lmax, mmax
#import TTCal: Nfreq

#include("special.jl")     # special functions
#include("physics.jl")     # physical constants and cosmology
#include("parallel.jl")    # tools for parallel processing
#include("definitions.jl") # defines all the types
#include("visibilities.jl")
#include("mmodes.jl")
#include("transfermatrix.jl")
#include("alm.jl")

##include("noise.jl")
#include("sky.jl")

end

