@testset "average-frequency-channels.jl" begin
    mmax  = 1
    frequencies = [45u"MHz", 74u"MHz", 100u"MHz", 150u"MHz"]
    bandwidth   = [24u"kHz", 24u"kHz",   1u"MHz",   2u"kHz"]

    input = MFBlockMatrix(NoFile(), mmax, frequencies, bandwidth)
    for β = 1:length(frequencies), m = 0:mmax
        input[m, β] = fill(complex(β, 1.0), 1, 1)
    end

    output = BPJSpec.average_frequency_channels(input, 2)
    @test BPJSpec.frequencies(output) == [59.5u"MHz", (100+150*0.002)/1.002*u"MHz"]
    @test BPJSpec.bandwidth(output)   == [  48u"kHz", 1.002u"MHz"]
    for m = 0:mmax
        @test output[m, 1] == fill(1.5+1.0im, 1, 1)
        @test output[m, 2] ≈ fill(3008/1002+1.0im, 1, 1)
    end
end

