struct TestBlockMatrix{S} <: BPJSpec.AbstractBlockMatrix{Vector{Int}, 1}
    storage :: S
    cache   :: BPJSpec.Cache{Vector{Int}}
    metadata1 :: Float64
    metadata2 :: String
end
BPJSpec.metadata_fields(matrix::TestBlockMatrix) = (matrix.metadata1, matrix.metadata2)
BPJSpec.nblocks(::Type{<:TestBlockMatrix}, metadata1, metadata2) = length(metadata2)
BPJSpec.linear_index(::TestBlockMatrix, idx) = idx-100
BPJSpec.indices(matrix::TestBlockMatrix) = 101:length(matrix.metadata2)+100

@testset "abstract-block-matrix.jl" begin
    matrix = create(TestBlockMatrix, 1.0, "hello")
    for idx = 101:100+length("hello")
        matrix[idx] = [idx, idx+1]
        @test matrix[idx] == [idx, idx+1]
    end
    @test length(matrix.cache) == length("hello")

    path = tempname()
    matrix1 = create(TestBlockMatrix, MultipleFiles(path), 3.14, "hi")
    try
        for idx = 101:100+length("hi")
            matrix1[idx] = rand(Int, 10)
        end
        matrix2 = BPJSpec.load(path)
        @test typeof(matrix2) == TestBlockMatrix{MultipleFiles}
        @test matrix1.metadata1 == matrix2.metadata1
        @test matrix1.metadata2 == matrix2.metadata2
        for idx = 101:100+length("hi")
            @test matrix1[idx] == matrix2[idx]
        end
        cache!(matrix2)
        matrix2[101] = rand(Int, 10)
        flush!(matrix2)
        @test matrix1[101] == matrix2[101]

        matrix3 = similar(matrix1)
        @test typeof(matrix3) == TestBlockMatrix{NoFile}
        @test matrix1.metadata1 == matrix3.metadata1
        @test matrix1.metadata2 == matrix3.metadata2

    finally
        rm(path, recursive=true, force=true)
    end
end

