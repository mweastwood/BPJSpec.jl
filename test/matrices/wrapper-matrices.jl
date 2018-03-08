@testset "wrapper-matrices.jl" begin
    mmax = 2
    frequencies = [74.0u"MHz", 100.0u"MHz"]
    bandwidth = [24u"kHz", 1.0u"MHz"]
    path = tempname()

    @testset "MBlockMatrix" begin
        for storage in (NoFile(), SingleFile(path), MultipleFiles(path))
            matrix = MBlockMatrix(storage, mmax)
            X = rand(Complex128, 5, 5)
            Y = rand(Complex128, 3, 3)
            matrix[0] = X
            matrix[2] = Y
            @test matrix[0] == X
            @test matrix[2] == Y
            rm(path, force=true, recursive=true)
        end
    end

    @testset "FBlockMatrix" begin
        for storage in (NoFile(), SingleFile(path), MultipleFiles(path))
            matrix = FBlockMatrix(storage, frequencies, bandwidth)
            X = rand(Complex128, 5, 5)
            Y = rand(Complex128, 3, 3)
            matrix[1] = X
            matrix[2] = Y
            @test matrix[1] == X
            @test matrix[2] == Y
            rm(path, force=true, recursive=true)
        end
    end

    @testset "MFBlockMatrix" begin
        for storage in (NoFile(), SingleFile(path), MultipleFiles(path))
            matrix = MFBlockMatrix(storage, mmax, frequencies, bandwidth)
            X = rand(Complex128, 5, 5)
            Y = rand(Complex128, 3, 3)
            matrix[0, 1] = X
            matrix[2, 2] = Y
            @test matrix[0, 1] == X
            @test matrix[2, 2] == Y
            rm(path, force=true, recursive=true)
        end
    end

    @testset "MBlockVector" begin
        for storage in (NoFile(), SingleFile(path), MultipleFiles(path))
            vector = MBlockVector(storage, mmax)
            X = rand(Complex128, 5)
            Y = rand(Complex128, 3)
            vector[0] = X
            vector[2] = Y
            @test vector[0] == X
            @test vector[2] == Y
            rm(path, force=true, recursive=true)
        end
    end

    @testset "FBlockVector" begin
        for storage in (NoFile(), SingleFile(path), MultipleFiles(path))
            vector = FBlockVector(storage, frequencies, bandwidth)
            X = rand(Complex128, 5)
            Y = rand(Complex128, 3)
            vector[1] = X
            vector[2] = Y
            @test vector[1] == X
            @test vector[2] == Y
            rm(path, force=true, recursive=true)
        end
    end

    @testset "MFBlockVector" begin
        for storage in (NoFile(), SingleFile(path), MultipleFiles(path))
            vector = MFBlockVector(storage, mmax, frequencies, bandwidth)
            X = rand(Complex128, 5)
            Y = rand(Complex128, 3)
            vector[0, 1] = X
            vector[2, 2] = Y
            @test vector[0, 1] == X
            @test vector[2, 2] == Y
            rm(path, force=true, recursive=true)
        end
    end
end

