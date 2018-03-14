struct TTCalMetadata
    frequencies   :: Vector{typeof(1.0*u"Hz")}
    times         :: Vector{Epoch}
    positions     :: Vector{Position}
    phase_centers :: Vector{Direction}
end

@testset "metadata.jl" begin

    frequencies = [74u"MHz", 76u"MHz"]
    times = [Epoch(epoch"UTC", 57365.5u"d")]
    positions = [Position(pos"WGS84", 1u"m", 0u"°", 0u"°"),
                 Position(pos"WGS84", 1u"m", 0u"°", 0.001u"°")]
    phase_centers = [Direction(dir"AZEL", 0u"°", 90u"°")]

    frame = ReferenceFrame()
    set!(frame, mean(positions))
    set!(frame, times[1])
    position = measure(frame, mean(positions), pos"ITRF")
    phase_center = measure(frame, phase_centers[1], dir"ITRF")

    ttcal_metadata = TTCalMetadata(frequencies, times, positions, phase_centers)
    bpjspec_metadata = BPJSpec.from_ttcal(ttcal_metadata)

    @test bpjspec_metadata.frequencies == frequencies
    @test bpjspec_metadata.bandwidth == [24.0u"kHz", 24.0u"kHz"]
    @test bpjspec_metadata.position ≈ position
    @test bpjspec_metadata.phase_center == phase_center
    @test bpjspec_metadata.baselines[1] == Baseline(baseline"ITRF", 0, 0, 0)
    @test bpjspec_metadata.baselines[2] == Baseline(baseline"ITRF", 0.0009649423882365227,
                                                    0.0, -110.57429326938583)
    @test bpjspec_metadata.baselines[3] == Baseline(baseline"ITRF", 0, 0, 0)

end

