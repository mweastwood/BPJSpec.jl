let lmax = 20, mmax = 10
    alm = Alm(Complex128,lmax,mmax)
    for m = 0:mmax
        N = lmax-m+1
        BPJSpec.setblock!(alm,fill(m,N),m)
        @test BPJSpec.getblock(alm,m) == fill(m,N)
    end
    for m = 0:mmax, l = m:lmax
        @test alm[l,m] == m
    end
end

let lmax = 20, mmax = 10, Nbase = 100, ν = 45e6
    B = TransferMatrix(Nbase,lmax,mmax,ν)
    for m = 0:mmax
        rand!(B[m].block)
    end

    alm = Alm(Complex128,lmax,mmax)
    rand!(alm.alm)
    v = B*alm

    alm′ = tikhonov(B,v,tolerance=0.0)
    @test alm.alm ≈ alm′.alm

    alm′ = tikhonov(B,v,tolerance=1e16)
    @test vecnorm(alm′.alm) < 1e-10
end

