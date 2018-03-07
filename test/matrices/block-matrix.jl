@testset "block-matrix.jl" begin
    @testset "no file" begin
        metadata = BPJSpec.MMax(2)
        matrix = BPJSpec.BlockMatrix{Matrix{Float64}, 1}(metadata)
        @test BPJSpec.used(matrix.cache) == true

        X = randn(5, 5)
        Y = randn(4, 4)
        matrix[3] = X
        @test matrix[3] == X
        matrix[3] = Y
        @test matrix[3] == Y
    end

    @testset "single file" begin
        path = tempname()
        metadata = BPJSpec.MMax(2)
        matrix = BPJSpec.BlockMatrix{Matrix{Float64}, 1}(SingleFile(path), metadata)
        @test BPJSpec.used(matrix.cache) == false

        X = randn(5, 5)
        Y = randn(4, 4)
        Z = randn(3, 3)
        matrix[1] = X
        matrix[2] = Y
        matrix[3] = Z
        @test matrix[1] == X
        @test matrix[2] == Y
        @test matrix[3] == Z
        cache!(matrix)
        @test matrix[1] == X
        @test matrix[2] == Y
        @test matrix[3] == Z
        matrix[1] = Y
        @test matrix[1] == Y

        matrix′ = BPJSpec.BlockMatrix{Matrix{Float64}, 1}(path)
        @test matrix′[1] == X
        @test matrix′[2] == Y
        @test matrix′[3] == Z

        rm(path, recursive=true)
    end

    @testset "multiple files" begin
        path = tempname()
        metadata = BPJSpec.MMax(2)
        matrix = BPJSpec.BlockMatrix{Matrix{Float64}, 1}(MultipleFiles(path), metadata)
        @test BPJSpec.used(matrix.cache) == false

        X = randn(5, 5)
        Y = randn(4, 4)
        Z = randn(3, 3)
        matrix[1] = X
        matrix[2] = Y
        matrix[3] = Z
        @test matrix[1] == X
        @test matrix[2] == Y
        @test matrix[3] == Z
        matrix[3] = Y
        @test matrix[3] == Y
        cache!(matrix)
        @test matrix[1] == X
        @test matrix[2] == Y
        @test matrix[3] == Y
        matrix[1] = Y
        @test matrix[1] == Y

        matrix′ = BPJSpec.BlockMatrix{Matrix{Float64}, 1}(path)
        @test matrix′[1] == X
        @test matrix′[2] == Y
        @test matrix′[3] == Y
        flush!(matrix)
        @test matrix′[1] == Y
        @test matrix′[2] == Y
        @test matrix′[3] == Y

        rm(path, recursive=true)
    end
end

