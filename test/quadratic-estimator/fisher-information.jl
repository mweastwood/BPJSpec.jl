@testset "fisher-information.jl" begin
    lmax = mmax = 2
    frequencies = [74.0u"MHz", 75.0u"MHz", 76.0u"MHz"]
    bandwidth   = [ 1.0u"MHz",  1.0u"MHz",  1.0u"MHz"]
    Nfreq = length(frequencies)
    Nbase = 10

    transfermatrix   = create(MBlockMatrix, mmax)
    covariancematrix = create(MBlockMatrix, mmax)
    basis = [create(LBlockMatrix, lmax, frequencies, bandwidth) for a = 1:2]

    for m = 0:mmax
        X = BPJSpec.two(m)*Nbase*Nfreq
        Y = Nfreq*(lmax - m + 1)
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

    iterations = get(ENV, "CI", false) ? 10 : 1000
    F = fisher_information(transfermatrix, covariancematrix, basis, iterations=iterations)

    # manually compute tr(C⁻¹ Ca C⁻¹ Cb)
    F′ = similar(F)
    lhs = create(MBlockMatrix, mmax)
    rhs = create(MBlockMatrix, mmax)
    for a = 1:length(basis), b = 1:length(basis)
        Ca = basis[a]
        Cb = basis[b]
        @. lhs = covariancematrix \ (transfermatrix * Ca * BPJSpec.T(transfermatrix))
        @. rhs = covariancematrix \ (transfermatrix * Cb * BPJSpec.T(transfermatrix))
        output = 0.0
        for m = 0:mmax
            output += real(sum(lhs[m] .* transpose(rhs[m])))
        end
        F′[a, b] = output
    end

    @test all(abs.(F .- F′) ./ abs.(F′) .< 0.2 * sqrt(1000/iterations))
end

