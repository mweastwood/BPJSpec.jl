using BPJSpec
using Base.Test
using Unitful, UnitfulAstro
using CasaCore.Measures

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

    include("spherical-harmonics.jl")

    @testset "physics" begin
        include("physics/cosmology.jl")
        include("physics/recombination-lines.jl")
    end

    @testset "sky" begin
        include("sky/foregrounds.jl")
        include("sky/signal.jl")
        include("sky/noise.jl")
    end

    @testset "matrices" begin
        include("matrices/storage-mechanisms.jl")
        include("matrices/block-matrix.jl")
    #    include("matrices/block-diagonal-matrix.jl")
    #    include("matrices/spectral-block-diagonal-matrix.jl")
    #    include("matrices/angular-covariance-matrix.jl")
    #    include("matrices/noise-covariance-matrix.jl")
    #    include("matrices/transfer-matrix.jl")
    end

    #@testset "algorithms" begin
    #    include("algorithms/average-frequency-channels.jl")
    #    #include("algorithms/full-rank-compress.jl")
    #end

    #@testset "quadratic-estimator" begin
    #    include("quadratic-estimator/mixing-matrix.jl")
    #end
end

