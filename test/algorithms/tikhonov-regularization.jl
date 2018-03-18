@testset "tikhonov-regularizatiion.jl" begin

    lmax = mmax = 1
    frequencies = [45u"MHz", 74u"MHz"]
    bandwidth   = [24u"kHz", 24u"kHz"]
    Nbase = 10
    Nfreq = length(frequencies)

    transfermatrix = create(MFBlockMatrix, mmax, frequencies, bandwidth)
    mmodes         = create(MFBlockVector, mmax, frequencies, bandwidth)
    alm            = create(MFBlockVector, mmax, frequencies, bandwidth)
    for β = 1:Nfreq, m = 0:mmax
        x = BPJSpec.two(m)*Nbase
        y = lmax - m + 1
        transfermatrix[m, β] = complex.(randn(x, y), randn(x, y))
    end

    @testset "regular" begin
        for β = 1:Nfreq, m = 0:mmax
            y = lmax - m + 1
            alm[m, β] = complex.(randn(y), randn(y))
        end
        @. mmodes = transfermatrix * alm

        alm′ = tikhonov(transfermatrix, mmodes, regularization=0.0, mfs=false)
        for β = 1:Nfreq, m = 0:mmax
            @test alm[m, β] ≈ alm′[m, β]
        end
    end

    @testset "mfs" begin
        for m = 0:mmax
            y = lmax - m + 1
            block = complex.(randn(y), randn(y))
            for β = 1:Nfreq
                alm[m, β] = block
            end
        end
        @. mmodes = transfermatrix * alm

        alm′ = tikhonov(transfermatrix, mmodes, regularization=0.0, mfs=true)
        for β = 1:Nfreq, m = 0:mmax
            @test alm[m, β] ≈ alm′[m]
        end
    end

end

