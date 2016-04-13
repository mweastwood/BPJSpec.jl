@testset "transfermatrix.jl" begin
    let Nbase = 10, lmax = 10
        @test BPJSpec.initial_block_size(TransferMatrix, Nbase, lmax, 0) == (10, 11)
        @test BPJSpec.initial_block_size(TransferMatrix, Nbase, lmax, 1) == (20, 10)
        @test BPJSpec.initial_block_size(TransferMatrix, Nbase, lmax, 2) == (20,  9)
    end

    let Nbase = 3, lmax = 7, mmax = 5, ν = 45e6
        transfermatrix = TransferMatrix(tempname(), lmax, mmax, [ν])
        BPJSpec.initialize!(transfermatrix, Nbase)
        for m = 0:mmax
            block = BPJSpec.read_block(transfermatrix, m, ν)
            @test all(block .== 0)
            @test size(block) == BPJSpec.initial_block_size(TransferMatrix, Nbase, lmax, m)

            oldblock = copy(block)
            rand!(oldblock)
            BPJSpec.write_block(transfermatrix, oldblock, m, ν)
            newblock = BPJSpec.read_block(transfermatrix, m, ν)
            @test oldblock == newblock
        end
    end

    # Test the plane wave expansion
    let u = 2.0, v = -3.0, w = 2.0
        alm_real, alm_imag = BPJSpec.planewave(u,v,w,0,100,100)
        map_real = alm2map(alm_real,512)
        map_imag = alm2map(alm_imag,512)
        map_test_real = HealpixMap(Float64,512)
        map_test_imag = HealpixMap(Float64,512)
        for i = 1:nside2npix(512)
            vec = LibHealpix.pix2vec_ring(512,i)
            phase = 2π*(u*vec[1]+v*vec[2]+w*vec[3])
            map_test_real[i] = cos(phase)
            map_test_imag[i] = sin(phase)
        end
        @test pixels(map_real) ≈ pixels(map_test_real)
        @test pixels(map_imag) ≈ pixels(map_test_imag)
    end

    #=
    # test that the elements for an autocorrelation are correct
    let Nbase = 1, lmax = 3, mmax = 3
        path = tempname()
        beam = ones(12*512*512) |> HealpixMap
        meta = Metadata([TTCal.Antenna(observatory("OVRO_MMA"))], [TTCal.Baseline(1, 1)], [45e6],
                        Direction(dir"AZEL", 0degrees, 90degrees), Epoch(epoch"UTC", (2015-1858)*365*days),
                        ConstantBeam())
        transfermatrix = gentransfermatrix(path, meta, lmax, mmax)

        for α = 1:Nbase
            @test abs(B[1][α,1]) > 1
        end
        for l = 1:lmax, α = 1:Nbase
            @test abs(B[1][α,l+1]) < 1e-5 # likely limited by accuracy of spherical harmonic transform
        end
        for m = 1:mmax, l = m:lmax, α = 1:2Nbase
            @test abs(B[m+1][α,l-m+1]) < 1e-5 # likely limited by accuracy of spherical harmonic transform
        end
    end
    =#

    #=
    # if we've defined the transfer matrix correctly, we should be able to compute
    # visibilities that match those computed directly from TTCal
    let Nant = 5, Nfreq = 2, lmax = 100, mmax = 100
        name,ms = createms(Nant,Nfreq)
        beam = HealpixMap(ones(12*512*512))
        u,v,w = itrf_baselines(ms)
        phasecenter = itrf_phasecenter(ms)

        # let's begin with a single source on the north pole
        alm = Alm(Complex128,lmax,mmax)
        for m = 0:mmax, l = m:lmax
            alm[l,m] = BPJSpec.Y(l,m,0.0,0.0) |> conj
        end
        B = transfer(beam,u,v,w,ms.ν[1],phasecenter,lmax=lmax,mmax=mmax)
        v = B*alm
        vis = visibilities(v)
        expected = ones(Complex128,size(vis)) # visibilities should always be unity
        @test isapprox(vis,expected,atol=1e-1)

        # note that there are two things likely contributing to the rough
        # tolerance:
        # 1. point sources carry power to large l, which means that by
        #    truncating the spherical harmonic expansion at some lmax,
        #    we are missing some of the flux
        # 2. alm2map and map2alm have some error that we must live with
        #    until I wrap their iterative counterparts
    end

    # test transfermatrix i/o
    let Nbase = 100, lmax = 20, mmax = 20
        filename = tempname()*".jld"
        ν = 45e6

        B1 = BPJSpec.TransferMatrix(Nbase,lmax,mmax,ν)
        for m = 0:mmax
            rand!(B1[m+1].block)
        end
        BPJSpec.save(filename,B1)

        B2 = BPJSpec.load(filename,lmax,mmax,ν)
        @test B1 == B2

        # and make sure we can write multiple frequencies to the same file
        B3 = BPJSpec.TransferMatrix(Nbase,lmax,mmax,ν+1e6)
        BPJSpec.save(filename,B3)

        B4 = BPJSpec.load(filename,lmax,mmax,ν)
        B5 = BPJSpec.load(filename,lmax,mmax,ν+1e6)
        @test B1 == B4
        @test B3 == B5
    end

    # test that we can make the transfer matrix square while leaving the
    # singular values untouched
    let Nbase = 100, lmax = 20, mmax = 20
        B1 = BPJSpec.TransferMatrix(Nbase,lmax,mmax,45e6)
        for m = 0:mmax
            rand!(B1[m+1].block)
        end

        P  = preserve_singular_values(B1)
        @test typeof(P) == BPJSpec.Blocks{BPJSpec.MatrixBlock,BPJSpec.NoMetadata}
        B2 = P*B1
        @test typeof(B2) == BPJSpec.TransferMatrix
        @test BPJSpec.is_single_frequency(B2)
        @test !BPJSpec.is_single_m(B2)

        for m = 0:mmax
            x,y = size(B2[m+1])
            @test x == y

            U1,S1,V1 = svd(B1[m+1])
            U2,S2,V2 = svd(B2[m+1])
            @test S1 ≈ S2
        end
        @test B1.meta == B2.meta
    end
    =#
end

