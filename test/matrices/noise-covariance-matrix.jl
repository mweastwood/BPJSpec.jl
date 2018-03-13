@testset "noise-covariance-matrix.jl" begin
    lmax = mmax = 1
    frequencies = [45u"MHz", 74u"MHz"]
    bandwidth   = [ 1u"MHz", 16u"MHz"]

    mmodes = BPJSpec.create(MModes, mmax, frequencies, bandwidth)
    for β = 1:length(frequencies), m = 0:mmax
        mmodes[m, β] = ones(Complex128, BPJSpec.two(m))
    end

    model = BPJSpec.NoiseModel(100u"K", 12u"hr", 2, 2.4u"sr")
    N = BPJSpec.create(NoiseCovarianceMatrix, mmax, frequencies, bandwidth)
    compute!(N, model, mmodes)

    for β = 1:length(frequencies), m = 0:mmax
        block = N[m, β]
        σ² = ustrip(uconvert(u"Jy^2", model(m, frequencies[β], bandwidth[β])^2))
        @test block == Diagonal(σ² .* ones(BPJSpec.two(m)))
    end
end

