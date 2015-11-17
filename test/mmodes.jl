let Nbase = 10, lmax = 10
    @test BPJSpec.default_size(BPJSpec.MModeBlock,Nbase,0) == (10,)
    @test BPJSpec.default_size(BPJSpec.MModeBlock,Nbase,1) == (20,)
    @test BPJSpec.default_size(BPJSpec.MModeBlock,Nbase,2) == (20,)

    @test (BPJSpec.MModeBlock(Nbase,0,45e6) |> size
                == BPJSpec.default_size(BPJSpec.MModeBlock,Nbase,0))
    @test (BPJSpec.MModeBlock(Nbase,1,45e6) |> size
                == BPJSpec.default_size(BPJSpec.MModeBlock,Nbase,1))
    @test (BPJSpec.MModeBlock(Nbase,2,45e6) |> size
                == BPJSpec.default_size(BPJSpec.MModeBlock,Nbase,2))
end

let Nbase = 10, lmax = 10, mmax = 10
    blocks = [BPJSpec.MModeBlock(Nbase,m,45e6) for m = 0:mmax]
    v = MModes(blocks)
    @test typeof(v) == MModes{BPJSpec.one_ν}

    blocks = [BPJSpec.MModeBlock(Nbase,3,ν) for ν in linspace(45e6,50e6,5)]
    v = MModes(blocks)
    @test typeof(v) == MModes{BPJSpec.one_m}
end

# m-mode generation and visibility generation should be inverse operations
let Nbase = 100, mmax = 100
    data   = rand(Complex128,Nbase,2mmax+1)
    mmodes = MModes(data,45e6,mmax=mmax)
    data′  = visibilities(mmodes)
    @test data ≈ data′
end

# test m-mode i/o
let Nbase = 100, mmax = 20
    filename = tempname()*".jld"
    ν = 45e6

    v1 = MModes(Nbase,mmax,ν)
    for m = 0:mmax
        rand!(v1[m].block)
    end
    save_mmodes(filename,v1)

    v2 = load_mmodes(filename,ν)
    @test v1 == v2

    # and make sure we can write multiple frequencies to the same file
    v3 = MModes(Nbase,mmax,ν+1e6)
    save_mmodes(filename,v3)

    v4 = load_mmodes(filename,ν)
    v5 = load_mmodes(filename,ν+1e6)
    @test v1 == v4
    @test v3 == v5
end

