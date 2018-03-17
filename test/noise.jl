# test NoiseModel
let
    Tsys0 = 6420
    α = 2.56
    ν0 = 45e6
    Δν = 24e3
    τ_total = 24*3600.0
    τ_int   = 30.0
    noise = NoiseModel(Tsys0,α,ν0,Δν,τ_total,τ_int)
    @test noise(0,ν0) ≈ Tsys0^2/(τ_total*Δν)
    @test noise(0,2ν0) ≈ Tsys0^2/(τ_total*Δν) * 2^(-2α)
end

# test that we can compress the noise covariance matrix in the same way as
# the transfer matrix and m-modes
let Nbase = 10, lmax = 2, mmax = 2
    B = BPJSpec.TransferMatrix(Nbase,lmax,mmax,45e6)
    for m = 0:mmax
        rand!(B[m+1].block)
    end
    P = preserve_singular_values(B)

    N = BPJSpec.DiagonalNoiseMatrix(Nbase,mmax,45e6)
    for m = 0:mmax
        for α = 1:size(N[m+1])[1]
            N[m+1][α,α] = rand(Complex128)
        end
    end

    N′ = P*N*P'
    @test typeof(N′) == BPJSpec.NoiseMatrix
    @test full(N′) ≈ full(P)*full(N)*full(P)'
end

