@testset "m-modes.jl" begin

    frequencies = [74.0u"MHz", 100.0u"MHz"]
    bandwidth = [24u"kHz", 1.0u"MHz"]

    frame = ReferenceFrame()
    position  = measure(frame, observatory("OVRO_MMA"), pos"ITRF")
    baselines = [Baseline(baseline"ITRF", 0, 0, 0)]
    phase_center = Direction(position)
    metadata  = BPJSpec.Metadata(frequencies, bandwidth, position, baselines, phase_center)
    hierarchy = BPJSpec.Hierarchy(metadata)

    mmodes = create(MModes, NoFile(), metadata, hierarchy)

    ϕ = linspace(0, 2π, 6629)[1:6628]
    X = reshape(cis.(ϕ) .+ 1, 6628, 1)

    compute!(MModes, mmodes, X, 1)

    @test mmodes[0, 1] ≈ [1]
    @test mmodes[1, 1] ≈ [1, 0]
    @test norm(mmodes[2, 1]) < eps(Float64)

end

