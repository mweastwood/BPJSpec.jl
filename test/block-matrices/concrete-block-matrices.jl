function test_matrix(T, N, idx1, idx2, fields...)
    path1 = tempname()
    path2 = tempname()
    for S in (NoFile, SingleFile, MultipleFiles)
        try
            matrix1 = BPJSpec.create(T, S(path1), fields...)
            matrix2 = BPJSpec.create(T, S(path2), fields...)
            if N == 1
                X = rand(Complex128, 5)
                Y = rand(Complex128, 3)
            else
                X = rand(Complex128, 5, 5)
                Y = rand(Complex128, 3, 3)
            end
            matrix1[idx1...] = X
            matrix1[idx2...] = Y
            @test matrix1[idx1...] == X
            @test matrix1[idx2...] == Y
            @. matrix2 = 2*(matrix1+matrix1)
            @test matrix2[idx1...] == 2*(X+X)
            @test matrix2[idx2...] == 2*(Y+Y)
            if S != NoFile
                matrix3 = BPJSpec.load(path1)
                @test matrix3[idx1...] == X
                @test matrix3[idx2...] == Y
                cache!(matrix1)
                cache!(matrix2)
                @. matrix2 = 3*(matrix1+matrix1)
                @test matrix1[idx1...] == X
                @test matrix1[idx2...] == Y
                @test matrix2[idx1...] == 3*(X+X)
                @test matrix2[idx2...] == 3*(Y+Y)
            end
        finally
            rm(path1, force=true, recursive=true)
            rm(path2, force=true, recursive=true)
        end
    end
end

@testset "concrete-block-matrices.jl" begin
    length = 2
    mmax   = 1
    frequencies = [74.0u"MHz", 100.0u"MHz"]
    bandwidth   = [24u"kHz", 1.0u"MHz"]

    @testset "SimpleBlockArray" begin
        test_matrix(SimpleBlockVector, 1, (1,), (2,), length)
        test_matrix(SimpleBlockMatrix, 2, (1,), (2,), length)
    end

    @testset "MBlockArray" begin
        test_matrix(MBlockVector, 1, (0,), (1,), mmax)
        test_matrix(MBlockMatrix, 2, (0,), (1,), mmax)
    end

    @testset "FBlockArray" begin
        test_matrix(FBlockVector, 1, (1,), (2,), frequencies, bandwidth)
        test_matrix(FBlockMatrix, 2, (1,), (2,), frequencies, bandwidth)
    end

    @testset "MFBlockArray" begin
        test_matrix(MFBlockVector, 1, (0, 1), (0, 2), 0, frequencies, bandwidth)
        test_matrix(MFBlockMatrix, 2, (0, 1), (0, 2), 0, frequencies, bandwidth)
    end
end

