@testset "average-frequency-channels.jl" begin
    path1 = tempname()
    path2 = tempname()
    mmax  = 1
    frequencies = [45u"MHz", 74u"MHz", 100u"MHz", 150u"MHz"]
    bandwidth   = [24u"kHz", 24u"kHz",   1u"MHz",   2u"kHz"]

    input = TransferMatrix(path1, mmax, frequencies, bandwidth)
    try
        for β = 1:length(frequencies), m = 0:mmax
            input[m, β] = fill(complex(β, 1.0), 1, 1)
        end
        output = BPJSpec.average_frequency_channels(input, 2, output=path2)
        try
            @test output.frequencies == [59.5u"MHz",   125u"MHz"]
            @test output.bandwidth   == [  48u"kHz", 1.002u"MHz"]
            for m = 0:mmax
                @test output[m, 1] == fill(1.5+1.0im, 1, 1)
                @test output[m, 2] == fill(3.5+1.0im, 1, 1)
            end
        finally
            rm(path2, recursive=true)
        end
    finally
        rm(path1, recursive=true)
    end
end

