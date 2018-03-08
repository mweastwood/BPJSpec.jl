@testset "mixing-matrix.jl" begin
    A = randn(5, 5)
    F = A'*A

    M⁻¹ = BPJSpec.inverse_mixing_matrix(F, strategy=:unwindowed)
    W   = BPJSpec.window_functions(F, M⁻¹)
    Σ   = BPJSpec.windowed_covariance(F, M⁻¹)
    @test norm(W - I) < 1e-13
    @test all(diag(Σ) .≥ inv.(diag(F)))

    M⁻¹ = BPJSpec.inverse_mixing_matrix(F, strategy=:minvariance)
    W   = BPJSpec.window_functions(F, M⁻¹)
    @test norm(M⁻¹ - diagm(diag(M⁻¹))) < 1e-13
    @test all(sum(W, 2) .≈ 1)

    M⁻¹ = BPJSpec.inverse_mixing_matrix(F, strategy=:uncorrelated)
    W   = BPJSpec.window_functions(F, M⁻¹)
    Σ   = BPJSpec.windowed_covariance(F, M⁻¹)
    @test all(sum(W, 2) .≈ 1)
    @test norm(Σ - diagm(diag(Σ))) < 1e-13
end

