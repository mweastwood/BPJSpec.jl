@testset "mmodes.jl" begin
    let Nfreq = 2, Nant = 5, mmax = 5
        Ntime = 2mmax+1
        Nbase = (Nant*(Nant+1))÷2
        path = tempname()
        meta = metadata(Nant, Nfreq)
        visibilities = GriddedVisibilities(tempname(), meta, Ntime)
        for β = 1:Nfreq
            visibilities.data[β][:] = 1
            visibilities.weights[β][:] = 1
        end
        mmodes = MModes(path, visibilities, mmax)
        @test mmodes.path == path
        @test mmodes.mmax == mmax
        @test mmodes.frequencies == meta.channels
        # given constant visibilities only m=0 should be nonzero
        for β = 1:Nfreq
            block = mmodes[0,β]
            @test block == ones(Complex128, Nbase)
            for m = 1:mmax
                block = mmodes[m,β]
                @test norm(block) < eps(Float64)
            end
        end

        # opening the m-modes again should see all the changes to the data
        for block in mmodes.blocks
            rand!(block)
        end
        mmodes′ = MModes(path)
        @test mmodes.path == mmodes′.path
        @test mmodes.mmax == mmodes′.mmax
        @test mmodes.frequencies == mmodes′.frequencies
        for β = 1:Nfreq, m = 0:mmax
            @test mmodes[m,β] == mmodes′[m,β]
        end
    end

    #=
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
    =#
end

