function constant_beam(azimuth, elevation)
    1.0
end
const constant_beam_solid_angle = 4π

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

    transfermatrix = create(TransferMatrix, path, metadata, simple_beam)
    lmax = mmax = transfermatrix.mmax

    @test transfermatrix.storage.hierarchy.divisions == [0, 32]
    @test transfermatrix.storage.hierarchy.baselines == [[1, 2, 3, 4]]

    @testset "creating beam maps" begin
        map  = BPJSpec.create_beam_map(simple_beam, metadata, (2, 4))
        @test size(map) == (2, 4)
        rhat = BPJSpec.unit_vectors(map)
        for jdx = 1:size(map, 2), idx = 1:size(map, 1)
            r = rhat[idx, jdx]
            elevation = asin(dot(r, up))
            @test map[idx, jdx] == max(0, sin(elevation))
        end
    end

    @testset "plane wave" begin
        # the amplitude of the plane wave should be unity everywhere
        map  = BPJSpec.create_beam_map(simple_beam, metadata, (2, 4))
        rhat = BPJSpec.unit_vectors(map)
        real_fringe, imag_fringe = BPJSpec.plane_wave(rhat, metadata.baselines[2],
                                                      metadata.phase_center,
                                                      u"c" / metadata.frequencies[1])
        for jdx = 1:size(map, 2), idx = 1:size(map, 1)
            @test hypot(real_fringe[idx, jdx], imag_fringe[idx, jdx]) ≈ 1
        end

        # if we put the phase center exactly on a pixel, that pixel should have zero imaginary
        # component
        real_fringe, imag_fringe = BPJSpec.plane_wave(rhat, metadata.baselines[2], rhat[1, 1],
                                                      u"c" / metadata.frequencies[1])
        @test real_fringe[1, 1] == 1
        @test imag_fringe[1, 1] == 0
    end

    @testset "uniform sky / auto-correlations" begin
        # here we check to see that a uniform 1 K sky produces m-modes with the right amplitude on
        # an auto-correlation

        alm = create(MFBlockVector, mmax, frequencies, bandwidth)
        for β = 1:length(frequencies)
            for m = 0:mmax
                block = zeros(Complex128, lmax - m + 1)
                if m == 0
                    block[1] = sqrt(4π) # 1 K constant sky brightness
                end
                alm[m, β] = block
            end
        end

        mmodes = create(MFBlockVector, mmax, frequencies, bandwidth)
        @. mmodes = transfermatrix * alm
        for β = 1:length(frequencies)
            block = mmodes[0, β]
            expected = ustrip(uconvert(u"Jy", simple_beam_solid_angle
                                        * 2u"k*K"*frequencies[β]^2/(u"c")^2))
            # Note the tolerance here is set by the resolution of the map used for the spherical
            # harmonic transform of the beam during transfer matrix generation.
            @test abs(block[1] - expected) < 0.05
        end
    end


    @testset "point source / cross-correlations" begin
        # check that we can create visibilities from a point source at the north pole correctly (a
        # point source at the north pole only impacts m=0, so it is easy to verify)

        BPJSpec.rm_old_blocks!(transfermatrix.storage)
        compute!(TransferMatrix, transfermatrix, metadata, constant_beam)

        alm = create(MFBlockVector, mmax, frequencies, bandwidth)
        for β = 1:length(frequencies)
            λ = u"c" / frequencies[β]
            factor = ustrip(uconvert(u"K", λ^2*u"Jy"/(2*u"k")))
            for m = 0:mmax
                block = zeros(Complex128, lmax - m + 1)
                for l = m:lmax
                    block[l - m + 1] = conj(BPJSpec.Y(l, m, 0, 0)) * factor
                end
                alm[m, β] = block
            end
        end

        mmodes = create(MFBlockVector, mmax, frequencies, bandwidth)
        @. mmodes = transfermatrix * alm
        for β = 1:length(frequencies)
            block = mmodes[0, β]
            λ = u"c" / frequencies[β]
            for α = 1:length(baselines)
                b = baselines[α]
                r = Direction(dir"ITRF", 0, 0, 1) - phase_center
                phase = uconvert(u"rad", 2π*dot(b, r)/λ)
                @test block[α] ≈ cis(phase)
            end
            for m = 1:mmax
                @test all(mmodes[m, β] .== 0)
            end
        end
    end

    rm(path, recursive=true)
end

