let N = 3
    A =  rand(Complex128,N,N)
    B =  rand(Complex128,N,N)
    C =  rand(Complex128,N,N)
    O = zeros(Complex128,N,N)

    blocks = BPJSpec.MatrixBlock[]
    push!(blocks,BPJSpec.MatrixBlock(A))
    push!(blocks,BPJSpec.MatrixBlock(B))
    push!(blocks,BPJSpec.MatrixBlock(C))
    block_diagonal = BPJSpec.Blocks(blocks)

    @test size(block_diagonal) == (3N,3N)
    @test full(block_diagonal) == [A O O;
                                   O B O;
                                   O O C]
end

let
    A =  rand(Complex128,3)
    B =  rand(Complex128,4)
    C =  rand(Complex128,5)

    blocks = BPJSpec.VectorBlock[]
    push!(blocks,BPJSpec.VectorBlock(A))
    push!(blocks,BPJSpec.VectorBlock(B))
    push!(blocks,BPJSpec.VectorBlock(C))
    block_vector = BPJSpec.Blocks(blocks)

    @test length(block_vector) == 3+4+5
    @test full(block_vector) == [A; B; C]
end

