@testset "noise-covariance-matrix.jl" begin
    path = tempname()
    lmax = mmax = 1
    frequencies = [45u"MHz", 74u"MHz"]
    bandwidth   = [ 1u"MHz", 16u"MHz"]

    frame = ReferenceFrame()
    position  = measure(frame, observatory("OVRO_MMA"), pos"ITRF")
    baselines = [Baseline(baseline"ITRF", 0, 0, 0)]
    phase_center = Direction(position)
    metadata  = BPJSpec.Metadata(frequencies, bandwidth, position, baselines, phase_center)
    hierarchy = BPJSpec.Hierarchy(metadata)

    model = BPJSpec.NoiseModel(100u"K", 12u"hr", 2, 2.4u"sr")
    N = create(NoiseCovarianceMatrix, path, model, metadata, hierarchy)

    for β = 1:length(frequencies), m = 0:mmax
        block = N[m, β]
        σ² = ustrip(uconvert(u"Jy^2", model(m, frequencies[β], bandwidth[β])^2))
        @test block == Diagonal(σ² .* ones(BPJSpec.two(m)))
    end

    rm(path, recursive=true)
end

