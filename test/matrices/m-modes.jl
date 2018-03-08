@testset "m-modes.jl" begin

    mmax = 2
    frequencies = [74.0u"MHz", 100.0u"MHz"]
    bandwidth = [24u"kHz", 1.0u"MHz"]

    mmodes = MModes(NoFile(), mmax, frequencies, bandwidth)

    ϕ = linspace(0, 2π, 6629)[1:6628]
    X = reshape(cis.(ϕ) .+ 1, 6628, 1)

    compute!(mmodes, X, 1)

    @test mmodes[0, 1] ≈ [1]
    @test mmodes[1, 1] ≈ [1, 0]
    @test norm(mmodes[2, 1]) < eps(Float64)

end

