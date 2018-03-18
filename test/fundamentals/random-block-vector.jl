@testset "random-block-vector.jl" begin

    @testset "WhiteNoiseBlockVector" begin
        w = WhiteNoiseBlockVector()
        @test w[1]    === BPJSpec.WhiteNoise()
        @test w[1, 2] === BPJSpec.WhiteNoise()
    end

    @testset "RandomBlockVector" begin
        size = (5, 5)
        Ntrials = 10_000

        C = create(MBlockMatrix, 0)
        C0 = complex.(randn(size), randn(size))
        C0 = C0*C0'
        @test isposdef(C0)
        C[0] = C0

        n = RandomBlockVector(C)
        C0′ = zeros(Complex128, size)
        for idx = 1:Ntrials
            n0 = n[0]
            C0′ .+= n0*n0'
        end
        C0′ ./= Ntrials

        @test norm(C0-C0′) / norm(C0) < 0.05
    end

end

