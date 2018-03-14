using BPJSpec
using Base.Test
using Unitful, UnitfulAstro
using CasaCore.Measures

srand(123)
@testset "BPJSpec Tests" begin
    @testset "wrappers" begin
        include("wrappers/FastTransformsWrapper.jl")
        include("wrappers/CosmologyWrapper.jl")
    end

    @testset "utilities" begin
        include("utilities/misc.jl")
        include("utilities/parallel.jl")
        include("utilities/recombination-lines.jl")
    end

    @testset "sky" begin
        include("sky/foregrounds.jl")
        include("sky/signal.jl")
        include("sky/noise.jl")
    end

    @testset "matrices" begin
        include("matrices/storage-mechanisms.jl")
        include("matrices/concrete-block-matrices.jl")
        include("matrices/broadcasting.jl")
        include("matrices/transfer-matrix.jl")
        include("matrices/noise-covariance-matrix.jl")
        include("matrices/angular-covariance-matrix.jl")
        include("matrices/m-modes.jl")
    end

    @testset "algorithms" begin
        include("algorithms/average-frequency-channels.jl")
        #include("algorithms/full-rank-compress.jl")
    end

    @testset "quadratic-estimator" begin
        include("quadratic-estimator/mixing-matrix.jl")
    end
end

