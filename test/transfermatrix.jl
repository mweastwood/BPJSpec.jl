@testset "transfermatrix.jl" begin
    #=
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

