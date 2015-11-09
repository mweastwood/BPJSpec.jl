let Nbase = 10, lmax = 10
    @test BPJSpec.default_size(BPJSpec.TransferMatrixBlock,Nbase,lmax,0) == (10,11)
    @test BPJSpec.default_size(BPJSpec.TransferMatrixBlock,Nbase,lmax,1) == (20,10)
    @test BPJSpec.default_size(BPJSpec.TransferMatrixBlock,Nbase,lmax,2) == (20, 9)

    @test (BPJSpec.TransferMatrixBlock(Nbase,lmax,0,45e6) |> size
                == BPJSpec.default_size(BPJSpec.TransferMatrixBlock,Nbase,lmax,0))
    @test (BPJSpec.TransferMatrixBlock(Nbase,lmax,1,45e6) |> size
                == BPJSpec.default_size(BPJSpec.TransferMatrixBlock,Nbase,lmax,1))
    @test (BPJSpec.TransferMatrixBlock(Nbase,lmax,2,45e6) |> size
                == BPJSpec.default_size(BPJSpec.TransferMatrixBlock,Nbase,lmax,2))
end

let Nbase = 10, lmax = 10, mmax = 10
    blocks = [BPJSpec.TransferMatrixBlock(Nbase,lmax,m,45e6) for m = 0:mmax]
    B = TransferMatrix(blocks)
    @test typeof(B) == TransferMatrix{one_ν}

    blocks = [BPJSpec.TransferMatrixBlock(Nbase,lmax,3,ν) for ν in linspace(45e6,50e6,5)]
    B = TransferMatrix(blocks)
    @test typeof(B) == TransferMatrix{one_m}
end

# infinitesimally short baselines should only have nonzero elements when l=m=0
let Nbase = 1, lmax = 3, mmax = 3
    beam = ones(12*512*512) |> HealpixMap
    u = rand(Nbase)*1e-16
    v = rand(Nbase)*1e-16
    w = zeros(Nbase)
    ν = 45e6
    phasecenter = (0,0,1)
    B = gentransfer(beam,u,v,w,ν,phasecenter,lmax=lmax,mmax=mmax)

    for α = 1:Nbase
        @test abs(B[α,0,0]) > 1
    end
    for l = 1:lmax, α = 1:Nbase
        @test abs(B[α,l,0]) < 1e-5 # likeley limited by accuracy of spherical harmonic transform
    end
    for m = 1:mmax, l = m:lmax, α = 1:2Nbase
        @test abs(B[α,l,m]) < 1e-5 # likeley limited by accuracy of spherical harmonic transform
    end

end

