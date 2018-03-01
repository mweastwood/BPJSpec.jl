@testset "noise.jl" begin
    Tsys = 1000u"K"
    ν  = 100u"MHz"
    Δν = 24u"kHz"
    τ  = 13u"s"
    N  = 1000 # number of integrations
    λ  = u"c" / ν
    Ω  = 2π*u"sr"
    A  = λ^2 / Ω
    @test BPJSpec.standard_error(Tsys, ν, Δν, τ, N, Ω) ≈
            u"k"*Tsys / (A * sqrt(Δν * τ * N))

    model = BPJSpec.NoiseModel(Tsys, τ, N, Ω)
    @test model(0, ν, Δν) == BPJSpec.standard_error(Tsys, ν, Δν, τ, N, Ω)
    @test abs(model(6628, ν, Δν)) < 1e-4u"Jy" # due to time smearing
end

