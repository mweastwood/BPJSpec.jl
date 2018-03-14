@testset "baseline-hierarchy.jl" begin

    frequencies = [45u"MHz", 74u"MHz"]
    bandwidth   = [24u"kHz", 24u"kHz"]

    frame = ReferenceFrame()
    position = measure(frame, observatory("OVRO_MMA"), pos"ITRF")
    up    = Direction(position)
    north = Direction(dir"ITRF", 0, 0, 1)
    north = Direction(north - dot(north, up)*up)
    east  = Direction(cross(north, up))
    baselines = [Baseline(baseline"ITRF", 0, 0, 0),
                 Baseline(baseline"ITRF",   10north.x,  10north.y,  10north.z),
                 Baseline(baseline"ITRF",  100east.x,  100east.y,  100east.z),
                 Baseline(baseline"ITRF", 1000up.x,   1000up.y,   1000up.z)]
    phase_center = up
    metadata = BPJSpec.Metadata(frequencies, bandwidth, position, baselines, phase_center)


    hierarchy = BPJSpec.compute_baseline_hierarchy(metadata, 50)
    @test repr(hierarchy) == """
    |--------------+-----------+-------------|
    | lmin to lmax | baselines |  disk space |
    |--------------+-----------+-------------|
    |    0 to   32 |         2 |    0.000 GB |
    |   32 to   50 |         0 |    0.000 GB |
    |--------------+-----------+-------------|
    |        total |         2 |    0.000 GB |
    |--------------+-----------+-------------|
    """
    @test BPJSpec.Nbase(hierarchy,  0) == 2
    @test BPJSpec.Nbase(hierarchy, 50) == 0

    hierarchy = BPJSpec.compute_baseline_hierarchy(metadata, 500)
    @test repr(hierarchy) == """
    |--------------+-----------+-------------|
    | lmin to lmax | baselines |  disk space |
    |--------------+-----------+-------------|
    |    0 to   32 |         2 |    0.000 GB |
    |   32 to  171 |         0 |    0.000 GB |
    |  171 to  500 |         1 |    0.007 GB |
    |--------------+-----------+-------------|
    |        total |         3 |    0.008 GB |
    |--------------+-----------+-------------|
    """
    @test BPJSpec.Nbase(hierarchy,   0) == 3
    @test BPJSpec.Nbase(hierarchy,  50) == 1
    @test BPJSpec.Nbase(hierarchy, 500) == 1

    hierarchy = BPJSpec.compute_baseline_hierarchy(metadata, 5000)
    @test repr(hierarchy) == """
    |--------------+-----------+-------------|
    | lmin to lmax | baselines |  disk space |
    |--------------+-----------+-------------|
    |    0 to   32 |         2 |    0.000 GB |
    |   32 to  171 |         0 |    0.000 GB |
    |  171 to 1207 |         1 |    0.043 GB |
    | 1207 to 1707 |         0 |    0.000 GB |
    | 1707 to 5000 |         1 |    0.745 GB |
    |--------------+-----------+-------------|
    |        total |         4 |    0.789 GB |
    |--------------+-----------+-------------|
    """
    @test BPJSpec.Nbase(hierarchy,    0) == 4
    @test BPJSpec.Nbase(hierarchy,   50) == 2
    @test BPJSpec.Nbase(hierarchy,  500) == 2
    @test BPJSpec.Nbase(hierarchy, 5000) == 1
end

