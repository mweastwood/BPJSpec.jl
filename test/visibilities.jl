@testset "visibilties.jl" begin
    let
        # test that the sidereal time is the same one day later
        sidereal_day = 86164.09054 # seconds
        meta = metadata(1, 1)
        meta_one_day_later = deepcopy(meta)
        meta_one_day_later.time = Epoch(epoch"UTC", meta.time.time+sidereal_day)
        # note the tolerance here is 0.01 seconds
        @test abs(BPJSpec.sidereal_time(meta) - BPJSpec.sidereal_time(meta_one_day_later)) < 1e-2/(24*3600)
    end

    let Nfreq = 2, Nant = 5, Ntime = 11
        Nbase = (Nant*(Nant+1))÷2
        path = tempname()
        meta = metadata(Nant, Nfreq)
        visibilities = GriddedVisibilities(path, meta, Ntime)
        @test visibilities.path == path
        @test visibilities.Nbase == Nbase
        @test visibilities.Ntime == Ntime
        @test visibilities.frequencies == meta.channels
        @test visibilities.origin == BPJSpec.sidereal_time(meta)
        @test length(visibilities.data) == Nfreq
        @test length(visibilities.weights) == Nfreq
        for β = 1:Nfreq
            @test all(visibilities.data[β] .== 0)
            @test all(visibilities.weights[β] .== 0)
            @test size(visibilities.data[β]) == (Nbase, Ntime)
            @test size(visibilities.weights[β]) == (Nbase, Ntime)
        end

        # opening the visibilities again should see any changes made
        # to the data and weights
        for β = 1:Nfreq
            rand!(visibilities.data[β])
            rand!(visibilities.weights[β])
        end
        visibilities′ = GriddedVisibilities(path)
        @test visibilities.path == visibilities′.path
        @test visibilities.Nbase == visibilities′.Nbase
        @test visibilities.Ntime == visibilities′.Ntime
        @test visibilities.frequencies == visibilities′.frequencies
        @test visibilities.origin == visibilities′.origin
        @test visibilities.data == visibilities′.data
        @test visibilities.weights == visibilities′.weights
    end

    let Nfreq = 2, Nant = 5, Ntime = 11
        # test that we can grid visibilities
        Nbase = (Nant*(Nant+1))÷2
        path = tempname()
        meta = metadata(Nant, Nfreq)
        gridded_visibilities = GriddedVisibilities(path, meta, Ntime)
        data = fill(JonesMatrix(1, rand(), rand(), 1), Nbase, Nfreq)
        flags = fill(false, Nbase, Nfreq)
        ungridded_visibilities = Visibilities(data, flags)
        grid!(gridded_visibilities, meta, ungridded_visibilities)
        expected_data = zeros(Complex128, Nbase, Ntime)
        expected_weights = zeros(Float64, Nbase, Ntime)
        expected_data[:,1] = 1
        expected_weights[:,1] = 1
        for β = 1:Nfreq
            @test gridded_visibilities.data[β] == expected_data
            @test gridded_visibilities.weights[β] == expected_weights
        end
    end

    let Nfreq = 2, Nant = 5, Ntime = 11
        # test that the visibilities are weighted correctly when using `getindex`
        Nbase = (Nant*(Nant+1))÷2
        path = tempname()
        meta = metadata(Nant, Nfreq)
        visibilities = GriddedVisibilities(path, meta, Ntime)
        for β = 1:Nfreq
            data = visibilities.data[β]
            weights = visibilities.weights[β]
            for time = 1:Ntime, α = 1:Nbase
                data[α,time] = complex(α + time, α + time)
                weights[α,time] = α + time
            end
        end
        for β = 1:Nfreq
            grid = visibilities[β]
            @test grid == fill(complex(1, 1), Nbase, Ntime)
        end
    end
end

