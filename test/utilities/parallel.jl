@testset "parallel.jl" begin
    # This is hard to test reliably on Travis, because you only get a single node to work with.
    # Therefore you can't test any of the multi-node parallelism we rely on. We'll just fake it here
    # to test some of the utilities though.

    hostname = chomp(readstring(`hostname`))
    workers  = BPJSpec.categorize_workers()
    @test length(workers.dict) == 1
    @test haskey(workers.dict, hostname)
    @test workers.dict[hostname] == [1]

    workers = BPJSpec.Workers(Dict("node1" => [1, 2, 3],
                                   "node2" => [4, 5],
                                   "node3" => [6]))
    @test repr(workers) == """
    | Workers
    |---------
    |   node1 : 1, 2, 3
    |   node2 : 4, 5
    |   node3 : 6
    """

    leaders = BPJSpec.leaders(workers)
    @test length(leaders) == 3
    @test 1 in leaders
    @test 4 in leaders
    @test 6 in leaders
end

