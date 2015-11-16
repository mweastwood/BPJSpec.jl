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

