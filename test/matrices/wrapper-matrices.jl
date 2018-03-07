@testset "wrapper-matrices.jl" begin
    mmax = 2
    frequencies = [74.0u"MHz", 100.0u"MHz"]
    bandwidth = [24u"kHz", 1.0u"MHz"]
    metadata1 = BPJSpec.MMax(mmax)
    metadata2 = BPJSpec.MMaxFrequencies(mmax, frequencies, bandwidth)

    @testset "MBlockMatrix" begin
        matrix = MBlockMatrix(NoFile(), metadata1)
        X = rand(Complex128, 5, 5)
        Y = rand(Complex128, 3, 3)
        matrix[0] = X
        matrix[2] = Y
        @test matrix[0] == X
        @test matrix[2] == Y
    end

    @testset "MFBlockMatrix" begin
        path = tempname()
        matrix = MFBlockMatrix(MultipleFiles(path), metadata2)
        X = rand(Complex128, 5, 5)
        Y = rand(Complex128, 3, 3)
        matrix[0, 1] = X
        matrix[2, 2] = Y
        @test matrix[0, 1] == X
        @test matrix[2, 2] == Y
        rm(path, force=true, recursive=true)
    end
end

