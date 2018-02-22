@testset "spectral-block-diagonal-matrix.jl" begin
    path1 = tempname()
    path2 = tempname()
    mmax  = 5
    frequencies = [45u"MHz", 74u"MHz"]
    bandwidth   = [ 1u"MHz",  1u"MHz"]
    A = DenseSpectralBlockDiagonalMatrix(path1, mmax, frequencies, bandwidth)
    B = DenseSpectralBlockDiagonalMatrix(path2, mmax, frequencies, bandwidth)
    try
        for β = 1:length(frequencies), m = 0:mmax
            A[m, β] = complex.(randn(5, 5), randn(5, 5))
        end

        @. B = A * A
        for β = 1:length(frequencies), m = 0:mmax
            Am = A[m, β]
            Bm = B[m, β]
            @test Bm == Am * Am
        end

        BPJSpec.cache!(A)
        BPJSpec.cache!(B)
        for β = 1:length(frequencies), m = 0:mmax
            Am = A[m, β]
            Bm = B[m, β]
            @test Bm == Am * Am
        end

        @. B = A + A
        for β = 1:length(frequencies), m = 0:mmax
            Am = A[m, β]
            Bm = B[m, β]
            @test Bm == Am + Am
        end

    finally
        rm(path1, recursive=true)
        rm(path2, recursive=true)
    end
end

