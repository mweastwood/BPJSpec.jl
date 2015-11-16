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

#=

# Test ForegroundModel
let
    ν0 = 45.0
    A = 123
    α = 1.0
    β = 1.0
    ζ = 1.0
    foreground = ForegroundModel(ν0,A,α,β,ζ)
    @test foreground(0,ν0,ν0) == A
end

# Test SphericalSignalModel
let
    k = logspace(log10(0.03),log10(0.3))
    P = ones(length(k)-1)
    signal = SphericalSignalModel(k,P)

    ν = 45.0
    z = BPJSpec.redshift(ν)
    r = BPJSpec.comoving_distance(z)
    @test BPJSpec.Csignal_spherical(0,ν,ν,k,P) ≈ (k[end]-k[1])/(π*r^2)
    @test signal(0,ν,ν) ≈ (k[end]-k[1])/(π*r^2)

    ν1 = 45.0
    ν2 = 45.1
    z1 = BPJSpec.redshift(ν1)
    z2 = BPJSpec.redshift(ν2)
    r1 = BPJSpec.comoving_distance(z1)
    r2 = BPJSpec.comoving_distance(z2)
    Δr = r2-r1
    @test(BPJSpec.Csignal_spherical(0,ν1,ν2,k,P)
            ≈ (sin(k[end]*Δr)-sin(k[1]*Δr))/(π*Δr*r1*r2))
    @test signal(0,ν1,ν2) ≈ (sin(k[end]*Δr)-sin(k[1]*Δr))/(π*Δr*r1*r2)
end

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

