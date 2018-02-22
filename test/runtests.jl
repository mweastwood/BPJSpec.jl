using BPJSpec
using Base.Test
using Unitful, UnitfulAstro
#using CasaCore.Measures
#using CasaCore.Tables
#using LibHealpix
#using TTCal
#using JLD

#include("setup.jl")

srand(123)
@testset "BPJSpec Tests" begin
    @testset "general" begin
        A = complex.(randn(5, 5), randn(5, 5))
        B = A*A'
        C = BPJSpec.fix(diagm([1.0, 0.2, -1e-16]))
        D = BPJSpec.fix([1.0 0.1; 0 -1.0])
        @test BPJSpec.T(A) == A'
        @test BPJSpec.T(B) == B
        @test BPJSpec.H(A) == 0.5*(A+A')
        @test BPJSpec.H(B) == B
        @test C == C'
        @test isposdef(C)
        @test D == D'
        @test !isposdef(D) # make sure we can't force arbitrary matrices to be positive definite
        @test BPJSpec.two( 0) == 1
        @test BPJSpec.two( 1) == 2
        @test BPJSpec.two(-1) == 2
    end

    include("cosmology.jl")
    include("spherical-harmonics.jl")

    @testset "sky" begin
        include("sky/foregrounds.jl")
        include("sky/signal.jl")
    end



end

