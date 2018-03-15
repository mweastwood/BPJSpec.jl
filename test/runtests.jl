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
        include("sky-components/foregrounds.jl")
        include("sky-components/signal.jl")
    end

    @testset "interferometer" begin
        include("interferometer/metadata.jl")
        include("interferometer/baseline-hierarchy.jl")
        include("interferometer/noise-model.jl")
    end

    @testset "block-matrices" begin
        include("block-matrices/storage-mechanisms.jl")
        include("block-matrices/abstract-block-matrix.jl")
        include("block-matrices/concrete-block-matrices.jl")
        include("block-matrices/random-block-vector.jl")
        include("block-matrices/broadcasting.jl")
    end

    @testset "fundamentals" begin
        include("fundamentals/transfer-matrix.jl")
        include("fundamentals/noise-covariance-matrix.jl")
        include("fundamentals/angular-covariance-matrix.jl")
        include("fundamentals/m-modes.jl")
    end

    @testset "algorithms" begin
        include("algorithms/average-frequency-channels.jl")
        #include("algorithms/full-rank-compress.jl")
    end

    @testset "quadratic-estimator" begin
        include("quadratic-estimator/mixing-matrix.jl")
    end
end

