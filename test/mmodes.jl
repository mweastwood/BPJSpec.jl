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

