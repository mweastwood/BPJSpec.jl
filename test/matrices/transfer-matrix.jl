@testset "transfer-matrix.jl" begin
    path = tempname()
    lmax = 2
    frequencies = [45u"MHz", 74u"MHz"]
    bandwidth   = [24u"kHz", 24u"kHz"]

    frame = ReferenceFrame()
    position = measure(frame, observatory("OVRO_MMA"), pos"ITRF")
    up    = Direction(position)
    north = Direction(dir"ITRF", 0, 0, 1)
    north = Direction(north - dot(north, up)*up)
    east  = Direction(cross(north, up))
    # one auto-correlation and one meter baselines in the three cardinal directions
    baselines = [Baseline(baseline"ITRF", 0, 0, 0),
                 Baseline(baseline"ITRF", north.x, north.y, north.z),
                 Baseline(baseline"ITRF",  east.x,  east.y,  east.z),
                 Baseline(baseline"ITRF",    up.x,    up.y,    up.z)]
    phase_center = up
    metadata = BPJSpec.Metadata(frequencies, bandwidth, position, baselines, phase_center)

    B = HierarchicalTransferMatrix(path, metadata)
    try
        @test B.hierarchy.divisions == [0, 2, 3]
        @test B.hierarchy.baselines == [[1], [2, 3, 4]]
    finally
        rm(path, recursive=true)
    end
end

