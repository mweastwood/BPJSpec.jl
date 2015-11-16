# test create_empty_visibilities
let Nfreq = 2, Nbase = 10, Ntime = 11
    frequencies = linspace(45e6,50e6,Nfreq)
    filename = tempname()*".h5"
    create_empty_visibilities(filename,Nbase,Ntime,frequencies)

    jldopen(filename,"r") do file
        @test file["Nfreq"] |> read == Nfreq
        @test file["Nbase"] |> read == Nbase
        @test file["Ntime"] |> read == Ntime
        for β = 1:Nfreq
            name  = @sprintf("%.3fMHz",frequencies[β]/1e6)
            @test name in names(file)
            group = file[name]
            @test group["data"] |> read == zeros(Complex128,Nbase,Ntime)
            @test group["weights"] |> read == zeros(Float64,Nbase,Ntime)
        end
    end
end

# test grid_visibilities
let Nfreq = 2, Nant = 5, Ntime = 11
    frequencies = linspace(45e6,50e6,Nfreq)
    Nbase = div(Nant*(Nant-1),2)
    filename = tempname()*".jld"
    create_empty_visibilities(filename,Nbase,Ntime,frequencies)

    data  = rand(Complex64,Nbase,Nfreq)
    flags = zeros(Bool,Nbase,Nfreq)
    grid_visibilities(filename,data,flags,frequencies,0.0)

    jldopen(filename,"r") do file
        for β = 1:Nfreq
            name  = @sprintf("%.3fMHz",frequencies[β]/1e6)
            @test name in names(file)
            group = file[name]
            @test squeeze(group["data"][:,1],2) == data[:,β]
            @test squeeze(group["weights"][:,1],2) == ones(Nbase)
            @test group["data"][:,2:Ntime] == zeros(Nbase,Ntime-1)
            @test group["weights"][:,2:Ntime] == zeros(Nbase,Ntime-1)
        end
    end

    for β = 1:Nfreq
        @test load_visibilities(filename,frequencies[β])[:,1] == data[:,β]
    end
end

