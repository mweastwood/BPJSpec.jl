@testset "karhunen-loeve-transforms.jl" begin
    lmax = mmax = 2
    frequencies = [74.0u"MHz", 75.0u"MHz", 76.0u"MHz"]
    bandwidth   = [ 1.0u"MHz",  1.0u"MHz",  1.0u"MHz"]
    Nfreq = length(frequencies)
    Nbase = 10

    alm                    = create(MBlockVector, mmax)
    input_mmodes           = create(MBlockVector, mmax)
    input_transfermatrix   = create(MBlockMatrix, mmax)
    input_noisematrix      = create(MBlockMatrix, mmax)
    input_signalmatrix     = create(LBlockMatrix, lmax, frequencies, bandwidth)
    input_foregroundmatrix = create(LBlockMatrix, lmax, frequencies, bandwidth)

    output_mmodes            = create(MBlockVector, mmax)
    output_transfermatrix    = create(MBlockMatrix, mmax)
    output_covariance        = create(MBlockMatrix, mmax)
    output_foreground_filter = create(MBlockMatrix, mmax)
    output_noise_whitener    = create(MBlockMatrix, mmax)

    for m = 0:mmax
        X = BPJSpec.two(m)*Nbase*Nfreq
        Y = Nfreq*(lmax - m + 1)
        input_transfermatrix[m] = complex.(randn(X, Y), randn(X, Y))

        N = complex.(randn(X, X), randn(X, X))
        N = N'*N
        input_noisematrix[m] = N

        alm[m] = complex.(randn(Y), randn(Y))
    end
    for l = 0:lmax
        C = randn(Nfreq, Nfreq)
        C = C'*C
        input_signalmatrix[L(l)] = C

        C = randn(Nfreq, Nfreq)
        C = C'*C
        input_foregroundmatrix[L(l)] = C
    end

    @. input_mmodes = input_transfermatrix * alm

    foreground_filter!(output_mmodes, output_transfermatrix, output_covariance,
                       output_foreground_filter, output_noise_whitener,
                       input_mmodes, input_transfermatrix, input_noisematrix,
                       input_signalmatrix, input_foregroundmatrix, threshold=1.0)

    B = input_transfermatrix
    V = output_foreground_filter
    W = output_noise_whitener

    # test that the 21-cm signal power is greater than the foreground power after filtering
    F = create(MBlockMatrix, mmax)
    S = create(MBlockMatrix, mmax)
    @. F = BPJSpec.fix(T(V) * B * input_foregroundmatrix * T(B) * V)
    @. S = BPJSpec.fix(T(V) * B * input_signalmatrix     * T(B) * V)
    for m = 0:mmax
        @test all(eigvals(F[m], S[m]) .< 1.0)
    end

    # test that the covariance matrix is diagonal
    for m = 0:mmax
        C = output_covariance[m]
        @test norm(C - diagm(diag(C))) < 1e-13 * norm(C)
    end

    # test that the m-modes were processed in the same way as the transfer matrix
    expected_mmodes = create(MBlockVector, mmax)
    @. expected_mmodes = output_transfermatrix * alm
    for m = 0:mmax
        @test output_mmodes[m] ≈ expected_mmodes[m]
    end
end

