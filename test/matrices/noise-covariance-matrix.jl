@testset "noise-covariance-matrix.jl" begin
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
end

