@testset "GSLWrapper.jl" begin
    for i = 1:10
        θ =  π*rand()
        ϕ = 2π*rand()
        @test BPJSpec.Y(0,  0, θ, ϕ) ≈ 1/2 * sqrt(1/π)
        @test BPJSpec.Y(1,  0, θ, ϕ) ≈ 1/2 * sqrt(3/π) * cos(θ)
        @test BPJSpec.Y(1, +1, θ, ϕ) ≈ -1/2 * sqrt(3/(2π)) * sin(θ) * cis(+ϕ)
        @test BPJSpec.Y(1, -1, θ, ϕ) ≈ +1/2 * sqrt(3/(2π)) * sin(θ) * cis(-ϕ)
    end
end

