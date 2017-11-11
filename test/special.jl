@testset "special.jl" begin
    # Test the spherical harmonics
    for i = 1:10
        θ =  π*rand()
        ϕ = 2π*rand()
        @test BPJSpec.Y(0,0,θ,ϕ) ≈ 1/2 * sqrt(1/π)
        @test BPJSpec.Y(1,0,θ,ϕ) ≈ 1/2 * sqrt(3/π) * cos(θ)
        @test BPJSpec.Y(1,+1,θ,ϕ) ≈ -1/2 * sqrt(3/(2π)) * sin(θ) * exp(+1im*ϕ)
        @test BPJSpec.Y(1,-1,θ,ϕ) ≈ +1/2 * sqrt(3/(2π)) * sin(θ) * exp(-1im*ϕ)
    end

    # Test the spherical Bessel function
    for i = 1:10
        x = 100*randn() |> abs
        @test BPJSpec.j(0,x) ≈ sin(x)/x
        @test BPJSpec.j(1,x) ≈ sin(x)/x^2 - cos(x)/x
        @test BPJSpec.j(2,x) ≈ (3/x^2 - 1)*sin(x)/x - 3cos(x)/x^2
    end

    # Test the matrix trace
    for i = 1:10
        A = rand(Float64,50,50)
        B = rand(Float64,50,50)
        @test trace(A*B) ≈ BPJSpec.tr(A,B)
        A = rand(Complex128,50,50)
        B = rand(Complex128,50,50)
        @test trace(A*B) ≈ BPJSpec.tr(A,B)
    end

    # Test force_hermitian and force_posdef
    let
        A = rand(10,10)
        @test ishermitian(BPJSpec.force_hermitian(A))
        B = diagm(linspace(1,-1e-11,10))
        @test isposdef(BPJSpec.force_posdef(B))
    end
end

