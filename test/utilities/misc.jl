@testset "misc.jl" begin
    @testset "matrix functions" begin
        A = complex.(randn(5, 5), randn(5, 5))
        B = A*A'
        C = BPJSpec.fix(diagm([1.0, 0.2, -1e-16]))
        D = BPJSpec.fix([1.0 0.1; 0 -1.0])
        @test BPJSpec.T(A) == A'
        @test BPJSpec.T(B) == B
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
        @test BPJSpec.L(1) + BPJSpec.L(2) == BPJSpec.L(3)
        @test collect(BPJSpec.L(1):BPJSpec.L(3)) == [BPJSpec.L(1), BPJSpec.L(2), BPJSpec.L(3)]
    end

    @testset "type utilities" begin
        @test BPJSpec.discard_type_parameters(Array{Int, 2}) == Array
        @test BPJSpec.discard_type_parameters(Dict{Float64, String}) == Dict
        @test BPJSpec.discard_type_parameters(BPJSpec.AbstractBlockMatrix{Matrix{Float64}, 2}) ==
                BPJSpec.AbstractBlockMatrix
    end
end
