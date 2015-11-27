using BPJSpec
using Base.Test
using CasaCore.Measures
using CasaCore.Tables
using LibHealpix
using TTCal
using JLD

include("setup.jl")

srand(123)
include("special.jl")
include("physics.jl")
include("blocks.jl")
include("itrf.jl")
include("visibilities.jl")
include("transfermatrix.jl")
include("mmodes.jl")
include("alm.jl")
include("noise.jl")
include("sky.jl")

