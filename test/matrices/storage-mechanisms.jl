@testset "storage-mechanisms.jl" begin
    @testset "no file" begin
        storage = NoFile()
        @test repr(storage) == "<no file>"
        @test BPJSpec.distribute_write(storage) == false
        @test BPJSpec.distribute_read(storage) == false
    end

    @testset "single file" begin
        path = tempname()
        storage = SingleFile(path)
        @test repr(storage) == path
        @test BPJSpec.distribute_write(storage) == false
        @test BPJSpec.distribute_read(storage) == true

        metadata = randn()
        BPJSpec.write_metadata(storage, metadata)
        storage′, metadata′ = BPJSpec.read_metadata(path)
        @test storage.path == storage′.path
        @test metadata == metadata′

        X = randn(10, 10)
        Y = randn(11, 11)
        storage[5] = X
        @test storage[5] == X
        storage[6, 7] = Y
        @test storage[6, 7] == Y

        rm(path, force=true, recursive=true)
    end

    @testset "multiple files" begin
        path = tempname()
        storage = MultipleFiles(path)
        @test repr(storage) == path
        @test BPJSpec.distribute_write(storage) == true
        @test BPJSpec.distribute_read(storage) == true

        metadata = randn()
        BPJSpec.write_metadata(storage, metadata)
        storage′, metadata′ = BPJSpec.read_metadata(path)
        @test storage.path == storage′.path
        @test metadata == metadata′

        X = randn(10, 10)
        Y = randn(11, 11)
        storage[5] = X
        @test storage[5] == X
        storage[6, 7] = Y
        @test storage[6, 7] == Y

        rm(path, force=true, recursive=true)
    end

    @testset "cache" begin
        cache = BPJSpec.Cache{Vector{Float64}}(3)
        @test cache.used[] == false
        @test typeof(cache.cache) == Vector{Vector{Float64}}
        @test size(cache.cache) == (3,)
        BPJSpec.set!(cache)
        @test cache.used[] == true
        BPJSpec.unset!(cache)
        @test cache.used[] == false
        X = randn(5)
        cache[2] = X
        @test cache[2] == cache.cache[2] == X
    end
end

