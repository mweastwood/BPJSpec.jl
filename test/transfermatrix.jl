@testset "transfermatrix.jl" begin
    # test the beam map
    let
        obs = Position(pos"ITRF", 10meters, 0degrees, 90degrees) # observatory is on the north pole
        time = Epoch(epoch"UTC", (2015-1858)*365*24*3600*seconds)
        antennas = [TTCal.Antenna(obs)]
        baselines = [TTCal.Baseline(1, 1)]
        channels = [45e6]
        phase_center = Direction(dir"ITRF", 0degrees, 90degrees)
        beam = SineBeam(1.0)
        meta = Metadata(antennas, baselines, channels, phase_center, time, beam)
        map = BPJSpec.beam_map(meta, channels[1])
        expected_map = HealpixMap(Float64, 512)
        for idx = 1:length(expected_map)
            θ,ϕ = LibHealpix.pix2ang_ring(512, idx)
            if θ < π/2
                expected_map[idx] = cos(θ)
            end
        end
        @test isring(map)
        @test nside(map) == 512
        @test norm(pixels(map)-pixels(expected_map)) / norm(pixels(expected_map)) < eps(Float64)
    end

    # test the plane wave expansion
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

    let Nfreq = 2, Nant = 3, mmax = 5
        Nbase = (Nant*(Nant+1))÷2
        lmax = mmax
        path = tempname()
        meta = metadata(Nant, Nfreq)
        transfermatrix = TransferMatrix(path, meta, lmax, mmax)
        @test transfermatrix.path == path
        @test transfermatrix.lmax == lmax
        @test transfermatrix.mmax == mmax
        @test transfermatrix.frequencies == meta.channels
        # the baselines are randomly oriented, so we can't test the values in the
        # transfer matrix blocks, but we can at least make sure they're all the
        # correct size
        for β = 1:Nfreq, m = 0:mmax
            block = transfermatrix[m,β]
            @test size(block) == (BPJSpec.two(m)*Nbase, lmax-m+1)
        end
        # note we will indirectly test that the matrix is calculated correctly
        # by verifying calculations made with the transfer matrix

        # opening the transfer matrix again should see all changes to the blocks
        for β = 1:Nfreq, m = 0:mmax
            block = transfermatrix[m,β]
            rand!(block)
            transfermatrix[m,β] = block
        end
        transfermatrix′ = TransferMatrix(path)
        @test transfermatrix == transfermatrix′
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

