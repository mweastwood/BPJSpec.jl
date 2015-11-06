# test create_empty_visibilities
let Nfreq = 2, Nbase = 10, Ntime = 11
    filename = tempname()*".h5"
    create_empty_visibilities(filename,Nfreq,Nbase,Ntime)

    jldopen(filename,"r") do file
        @test file["Nfreq"] |> read == Nfreq
        @test file["Nbase"] |> read == Nbase
        @test file["Ntime"] |> read == Ntime
        for β = 1:Nfreq
            @test string(β) in names(file)
            @test file["$β/data"] |> read == zeros(Complex128,Nbase,Ntime)
            @test file["$β/weights"] |> read == zeros(Float64,Nbase,Ntime)
        end
    end
end

# test grid_visibilities
let Nfreq = 2, Nant = 5, Ntime = 11
    Nbase = div(Nant*(Nant-1),2)
    filename = tempname()*".jld"
    create_empty_visibilities(filename,Nfreq,Nbase,Ntime)

    data  = rand(Complex64,Nbase,Nfreq)
    flags = zeros(Bool,Nbase,Nfreq)
    grid_visibilities(filename,0.0,data,flags)

    jldopen(filename,"r") do file
        for β = 1:Nfreq
            @test squeeze(file["$β/data"][:,1],2) == data[:,β]
            @test squeeze(file["$β/weights"][:,1],2) == ones(Nbase)
            @test file["$β/data"][:,2:Ntime] == zeros(Nbase,Ntime-1)
            @test file["$β/weights"][:,2:Ntime] == zeros(Nbase,Ntime-1)
        end
    end

    for β = 1:Nfreq
        @test load_visibilities(filename,β)[:,1] == data[:,β]
    end
end

