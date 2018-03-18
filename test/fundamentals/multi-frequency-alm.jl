@testset "multi-frequency-alm.jl" begin
    lmax = mmax = 1
    frequencies = [74.0u"MHz", 100.0u"MHz"]
    bandwidth   = [24.0u"kHz",   1.0u"MHz"]

    input = create(MFBlockVector, mmax, frequencies, bandwidth)
    for β = 1:length(frequencies), m = 0:mmax
        block = zeros(Complex128, lmax - m + 1)
        for l = m:lmax
            block[l - m + 1] = complex(l*m, β)
        end
        input[m, β] = block
    end
    output = create(MultiFrequencyAlm, input)
    for m = 0:mmax, l = m:lmax
        block = zeros(Complex128, length(frequencies))
        for β = 1:length(frequencies)
            block[β] = complex(l*m, β)
        end
        @test output[l, m] == block
    end

    input = create(MBlockVector, mmax)
    for m = 0:mmax
        block = zeros(Complex128, length(frequencies)*(lmax - m + 1))
        idx = 1
        for β = 1:length(frequencies), l = m:lmax
            block[idx] = complex(l*m, β)
            idx += 1
        end
        input[m] = block
    end
    output = create(MultiFrequencyAlm, input, frequencies, bandwidth)
    for m = 0:mmax, l = m:lmax
        block = zeros(Complex128, length(frequencies))
        for β = 1:length(frequencies)
            block[β] = complex(l*m, β)
        end
        @test output[l, m] == block
    end

end

