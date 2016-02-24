let Nbase = 10
    @test BPJSpec.initial_block_size(BPJSpec.MModes,Nbase,0) == (10,)
    @test BPJSpec.initial_block_size(BPJSpec.MModes,Nbase,1) == (20,)
    @test BPJSpec.initial_block_size(BPJSpec.MModes,Nbase,2) == (20,)

    v = BPJSpec.MModes(Nbase,2,45e6)
    @test length(v[1]) == 10
    @test length(v[2]) == 20
    @test length(v[3]) == 20
end

# m-mode generation and visibility generation should be inverse operations
let Nbase = 100, mmax = 100
    data = rand(Complex128,Nbase,2mmax+1)
    v = mmodes(data;frequency=45e6,mmax=mmax)
    data′ = visibilities(v)
    @test data ≈ data′
end

# test m-mode i/o
let Nbase = 100, mmax = 20
    filename = tempname()*".jld"
    ν = 45e6

    v1 = BPJSpec.MModes(Nbase,mmax,ν)
    for m = 0:mmax
        rand!(v1[m+1].block)
    end
    BPJSpec.save(filename,v1)

    v2 = BPJSpec.load(filename,mmax,ν)
    @test v1 == v2

    # and make sure we can write multiple frequencies to the same file
    v3 = BPJSpec.MModes(Nbase,mmax,ν+1e6)
    BPJSpec.save(filename,v3)

    v4 = BPJSpec.load(filename,mmax,ν)
    v5 = BPJSpec.load(filename,mmax,ν+1e6)
    @test v1 == v4
    @test v3 == v5
end

# test that we can compress the m-modes without affecting the spherical
# harmonic coefficients
let Nbase = 100, lmax = 20, mmax = 20
    B1 = BPJSpec.TransferMatrix(Nbase,lmax,mmax,45e6)
    for m = 0:mmax
        rand!(B1[m+1].block)
    end

    v1 = BPJSpec.MModes(Nbase,mmax,45e6)
    for m = 0:mmax
        rand!(v1[m+1].block)
    end

    P  = preserve_singular_values(B1)
    B2 = P*B1
    v2 = P*v1

    @test typeof(v2) == BPJSpec.MModes
    @test BPJSpec.is_single_frequency(v2)
    @test !BPJSpec.is_single_m(v2)

    @test tikhonov(B1,v1,tolerance=0.0).alm ≈ tikhonov(B2,v2,tolerance=0.0).alm
    @test v1.meta == v2.meta
end

