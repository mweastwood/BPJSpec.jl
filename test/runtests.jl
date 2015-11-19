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
include("itrf.jl")
include("visibilities.jl")
include("transfermatrix.jl")
include("mmodes.jl")
include("alm.jl")
include("covariancematrix.jl")

#=
# Test SpectralMModes
let
    mmax = 10
    a = BPJSpec.MModeBlock{mmax}(0,complex(randn(5),randn(5)))
    b = BPJSpec.MModeBlock{mmax}(0,complex(randn(5),randn(5)))
    v = SpectralMModes{mmax}(0,[a,b])
    @test full(v) == [a.block; b.block]
end

# Test SpectralTransferMatrix
let
    lmax = 10
    mmax = 10
    a = BPJSpec.TransferMatrixBlock{lmax,mmax}(0,complex(randn(5,5),randn(5,5)))
    b = BPJSpec.TransferMatrixBlock{lmax,mmax}(0,complex(randn(5,5),randn(5,5)))
    z = zeros(Complex128,5,5)
    B = SpectralTransferMatrix{lmax,mmax}(0,[a,b])
    @test full(B) == [a.block z;
                      z b.block]
end
=#

