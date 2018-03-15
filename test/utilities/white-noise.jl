@testset "white-noise.jl" begin

    N = 1000
    w = BPJSpec.WhiteNoise()
    v = w(N)

    @test abs(v'*w)/N < 5/sqrt(N)
    @test abs(mean(v+w)) < 5/sqrt(N)
    @test abs(mean(w+v)) < 5/sqrt(N)
    @test abs(mean(v-w)) < 5/sqrt(N)
    @test abs(mean(w-v)) < 5/sqrt(N)

end

