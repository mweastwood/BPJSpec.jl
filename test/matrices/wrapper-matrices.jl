@testset "wrapper-matrices.jl" begin
    mmax = 2
    frequencies = [74.0u"MHz", 100.0u"MHz"]
    bandwidth = [24u"kHz", 1.0u"MHz"]

    @testset "MBlockMatrix" begin
        matrix = MBlockMatrix(NoFile(), mmax)
        X = rand(Complex128, 5, 5)
        Y = rand(Complex128, 3, 3)
        matrix[0] = X
        matrix[2] = Y
        @test matrix[0] == X
        @test matrix[2] == Y
    end

    @testset "MFBlockMatrix" begin
        path = tempname()
        matrix = MFBlockMatrix(MultipleFiles(path), mmax, frequencies, bandwidth)
        X = rand(Complex128, 5, 5)
        Y = rand(Complex128, 3, 3)
        matrix[0, 1] = X
        matrix[2, 2] = Y
        @test matrix[0, 1] == X
        @test matrix[2, 2] == Y
        rm(path, force=true, recursive=true)
    end
end

