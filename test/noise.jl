let Nbase = 10
    @test BPJSpec.initial_block_size(BPJSpec.DiagonalNoiseMatrix,Nbase,0) == (10,)
    @test BPJSpec.initial_block_size(BPJSpec.DiagonalNoiseMatrix,Nbase,1) == (20,)
    @test BPJSpec.initial_block_size(BPJSpec.DiagonalNoiseMatrix,Nbase,2) == (20,)

    N = BPJSpec.DiagonalNoiseMatrix(Nbase,2,45e6)
    @test size(N[1]) == (10,10)
    @test size(N[2]) == (20,20)
    @test size(N[3]) == (20,20)
end

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

# test noise covariance matrix i/o
let Nbase = 100, mmax = 20
    filename = tempname()*".jld"
    ν = 45e6

    N1 = BPJSpec.NoiseMatrix([BPJSpec.MatrixBlock(rand(Complex128,Nbase,Nbase)) for m = 0:mmax],
                             BPJSpec.NoiseMeta(mmax,ν))
    BPJSpec.save(filename,N1)

    N2 = BPJSpec.load(filename,mmax,ν)
    @test N1 == N2

    # and make sure we can write multiple frequencies to the same file
    N3 = BPJSpec.NoiseMatrix([BPJSpec.MatrixBlock(rand(Complex128,Nbase,Nbase)) for m = 0:mmax],
                             BPJSpec.NoiseMeta(mmax,ν+1e6))
    BPJSpec.save(filename,N3)

    N4 = BPJSpec.load(filename,mmax,ν)
    N5 = BPJSpec.load(filename,mmax,ν+1e6)
    @test N1 == N4
    @test N3 == N5
end

