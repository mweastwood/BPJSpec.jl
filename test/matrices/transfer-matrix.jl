function simple_beam(azimuth, elevation)
    if elevation > 0
        return sin(elevation)
    else
        return 0.0
    end
end
const simple_beam_solid_angle = π

@testset "transfer-matrix.jl" begin
    path = tempname()
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

    transfermatrix = TransferMatrix(path, metadata)
    try
        @test transfermatrix.storage.hierarchy.divisions == [0, 32]
        @test transfermatrix.storage.hierarchy.baselines == [[1, 2, 3, 4]]
        BPJSpec.compute!(transfermatrix, simple_beam)

        ############################################################################################
        # Construct a uniform sky and verify that the auto-correlation gets the correct amplitude.
        alm = MFBlockVector(NoFile(), transfermatrix.mmax, frequencies, bandwidth)
        for β = 1:length(frequencies)
            for m = 0:transfermatrix.mmax
                block = zeros(Complex128, transfermatrix.lmax - m + 1)
                if m == 0
                    block[1] = sqrt(4π) # 1 K constant sky brightness
                end
                alm[m, β] = block
            end
        end

        mmodes = MFBlockVector(NoFile(), transfermatrix.mmax, frequencies, bandwidth)
        @. mmodes = transfermatrix * alm
        for β = 1:length(frequencies)
            block = mmodes[0, β]
            expected = ustrip(uconvert(u"Jy", simple_beam_solid_angle
                                       * 2u"k*K"*frequencies[β]^2/(u"c")^2))
            # Note the tolerance here is set by the resolution of the map used for the spherical
            # harmonic transform of the beam during transfer matrix generation.
            @test abs(block[1] - expected) < 0.05
        end

    finally
        rm(path, recursive=true)
    end
end

