@testset "q-estimator.jl" begin
    lmax = mmax = 2
    frequencies = [74.0u"MHz", 75.0u"MHz", 76.0u"MHz"]
    bandwidth   = [ 1.0u"MHz",  1.0u"MHz",  1.0u"MHz"]
    Nfreq = length(frequencies)
    Nbase = 10

    mmodes           = create(MBlockVector, mmax)
    transfermatrix   = create(MBlockMatrix, mmax)
    covariancematrix = create(MBlockMatrix, mmax)
    basis = [create(LBlockMatrix, lmax, frequencies, bandwidth) for a = 1:2]

    for m = 0:mmax
        X = BPJSpec.two(m)*Nbase*Nfreq
        Y = Nfreq*(lmax - m + 1)
        mmodes[m]         = complex.(randn(X),    randn(X))
        transfermatrix[m] = complex.(randn(X, Y), randn(X, Y))

        N = complex.(randn(X, X), randn(X, X))
        N = N'*N
        covariancematrix[m] = N
    end

    for l = 0:lmax
        for a = 1:length(basis)
            C = randn(Nfreq, Nfreq)
            C = C'*C
            basis[a][BPJSpec.L(l)] = C
        end
    end

    q = q_estimator(mmodes, transfermatrix, covariancematrix, basis)
    for a = 1:length(basis)
        expected = 0.0
        for m = 0:mmax
            v = mmodes[m]
            B = transfermatrix[m]
            C = covariancematrix[m]
            Ca = B*basis[a][m]*B'
            w  = C\v
            expected += real(w'*Ca*w)
        end
        @test expected ≈ q[a]
    end

end

