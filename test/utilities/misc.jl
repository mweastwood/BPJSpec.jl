@testset "misc.jl" begin
    @testset "matrix functions" begin
        A = complex.(randn(5, 5), randn(5, 5))
        B = A*A'
        C = BPJSpec.fix(diagm([1.0, 0.2, -1e-16]))
        D = BPJSpec.fix([1.0 0.1; 0 -1.0])
        @test T(A) == A'
        @test T(B) == B
        @test BPJSpec.H(A) == 0.5*(A+A')
        @test BPJSpec.H(B) == B
        @test C == C'
        @test isposdef(C)
        @test D == D'
        @test !isposdef(D) # make sure we can't force arbitrary matrices to be positive definite
    end

    @testset "two" begin
        @test BPJSpec.two( 0) == 1
        @test BPJSpec.two( 1) == 2
        @test BPJSpec.two(-1) == 2
    end

    @testset "L" begin
        @test L(1) + L(2) == L(3)
        @test collect(L(1):L(3)) == [L(1), L(2), L(3)]
    end

    @testset "type utilities" begin
        @test BPJSpec.discard_type_parameters(Array{Int, 2}) == Array
        @test BPJSpec.discard_type_parameters(Dict{Float64, String}) == Dict
        @test BPJSpec.discard_type_parameters(BPJSpec.AbstractBlockMatrix{Matrix{Float64}, 2}) ==
                BPJSpec.AbstractBlockMatrix
    end

    @testset "stacking diagonally" begin
        A = complex.(randn(5, 6), randn(5, 6))
        B = complex.(randn(4, 3), randn(4, 3))
        Z1 = zeros(size(A, 1), size(B, 2))
        Z2 = zeros(size(B, 1), size(A, 2))
        @test BPJSpec.stack_diagonally([A, B]) == [A Z1; Z2 B]

        a = complex.(randn(5), randn(5))
        b = complex.(randn(4), randn(4))
        @test BPJSpec.stack_diagonally([a, b]) == [a; b]
    end
end

