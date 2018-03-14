@testset "foregrounds.jl" begin
    ν0 = 45u"MHz"
    A  = 123u"K^2"
    α  = 1.0
    β  = 1.0
    ζ  = Inf
    model = BPJSpec.ForegroundComponent(ν0, A, α, β, ζ)
    @test model(1000, ν0, ν0) == A
    @test model(20, ν0, ν0) / model(10, ν0, ν0) ≈ ((20+1)/(10+1))^(-α)
    @test model(0, 2ν0, 2ν0) / model(0, ν0, ν0) ≈ 2^(-2β)

    model = BPJSpec.extragalactic_point_sources()
    @test model(1000, 130u"MHz", 130u"MHz") == 57.0u"mK^2"
    model = BPJSpec.extragalactic_free_free()
    @test model(1000, 130u"MHz", 130u"MHz") == 0.014u"mK^2"
    model = BPJSpec.galactic_synchrotron()
    @test model(1000, 130u"MHz", 130u"MHz") == 700u"mK^2"
    model = BPJSpec.galactic_free_free()
    @test model(1000, 130u"MHz", 130u"MHz") == 0.088u"mK^2"

    str = "ForegroundComponent(ν0 = 130.000 MHz, A = 700.000 mK², α = 2.400, β = 2.800, ζ = 4.000)"
    @test repr(BPJSpec.galactic_synchrotron()) == str
end

