@testset "full-rank-compress.jl" begin

    lmax = mmax = 1
    frequencies = [45u"MHz", 74u"MHz"]
    bandwidth   = [24u"kHz", 24u"kHz"]
    Nbase = 100
    Nfreq = length(frequencies)

    transfermatrix = create(MFBlockMatrix, mmax, frequencies, bandwidth)
    mmodes         = create(MFBlockVector, mmax, frequencies, bandwidth)
    alm            = create(MFBlockVector, mmax, frequencies, bandwidth)
    noisematrix    = create(MFBlockMatrix, mmax, frequencies, bandwidth)
    for β = 1:Nfreq, m = 0:mmax
        x = BPJSpec.two(m)*Nbase
        y = lmax - m + 1
        transfermatrix[m, β] = complex.(randn(x, y), randn(x, y))
        alm[m, β] = complex.(randn(y), randn(y))
        block = complex.(randn(x, x), randn(x, x))
        noisematrix[m, β] = block'*block
    end
    @. mmodes = transfermatrix * alm

    transfermatrix′ = create(MFBlockMatrix, mmax, frequencies, bandwidth)
    mmodes′         = create(MFBlockVector, mmax, frequencies, bandwidth)
    noisematrix′    = create(MFBlockMatrix, mmax, frequencies, bandwidth)
    full_rank_compress!(mmodes′, transfermatrix′, noisematrix′,
                        mmodes,  transfermatrix,  noisematrix)

    # double check that the transfer matrix is the right size and that the singular values are
    # unchanged
    for β = 1:Nfreq, m = 0:mmax
        N = lmax - m + 1
        @test size(transfermatrix′[m, β]) == (N, N)
        U1, S1, V1 = svd(transfermatrix[m, β])
        U2, S2, V2 = svd(transfermatrix′[m, β])
        @test S1 ≈ S2
    end

    # double check that the transfer matrix and m-mode compressions are consistent
    mmodes″ = create(MFBlockVector, mmax, frequencies, bandwidth)
    @. mmodes″ = transfermatrix′ * alm
    for β = 1:Nfreq, m = 0:mmax
        @test mmodes′[m, β] ≈ mmodes″[m, β]
    end

    # double check that the noise matrix is still positive definite
    for β = 1:Nfreq, m = 0:mmax
        @test isposdef(noisematrix′[m, β])
    end
end

