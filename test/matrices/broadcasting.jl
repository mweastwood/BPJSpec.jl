@testset "broadcasting.jl" begin
    @testset "MBlockMatrix" begin
        path1 = tempname()
        path2 = tempname()
        mmax = 5
        A = MBlockMatrix(MultipleFiles(path1), mmax)
        B = MBlockMatrix(MultipleFiles(path2), mmax)
        try
            for m = 0:mmax
                A[m] = complex.(randn(5, 5), randn(5, 5))
            end

            @. B = A * A
            for m = 0:mmax
                Am = A[m]
                Bm = B[m]
                @test Bm == Am * Am
            end

            cache!(A)
            cache!(B)
            for m = 0:mmax
                Am = A[m]
                Bm = B[m]
                @test Bm == Am * Am
            end

            @. B = A + A
            for m = 0:mmax
                Am = A[m]
                Bm = B[m]
                @test Bm == Am + Am
            end

        finally
            rm(path1, force=true, recursive=true)
            rm(path2, force=true, recursive=true)
        end
    end

    @testset "MFBlockMatrix" begin
        path1 = tempname()
        path2 = tempname()
        mmax  = 5
        frequencies = [45u"MHz", 74u"MHz"]
        bandwidth   = [ 1u"MHz",  1u"MHz"]
        A = MFBlockMatrix(SingleFile(path1), mmax, frequencies, bandwidth)
        B = MFBlockMatrix(SingleFile(path2), mmax, frequencies, bandwidth)
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

            cache!(A)
            cache!(B)
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
end

