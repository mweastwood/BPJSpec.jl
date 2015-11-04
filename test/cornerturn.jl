# test touch
let Nfreq = 2, Nbase = 10, Ntime = 11
    filename = tempname()*".h5"
    BPJSpec.touch(filename,Nfreq,Nbase,Ntime)

    h5open(filename,"r") do file
        @test attrs(file)["Nfreq"] |> read == Nfreq
        @test attrs(file)["Nbase"] |> read == Nbase
        @test attrs(file)["Ntime"] |> read == Ntime
        for β = 1:Nfreq
            @test string(β) in names(file)
            group = file[string(β)]
            @test group["data"] |> read == zeros(2,Nbase,Ntime)
            @test group["weights"] |> read == zeros(Nbase,Ntime)
        end
    end
end

# test grid
let Nfreq = 2, Nant = 5, Ntime = 11
    Nbase = div(Nant*(Nant-1),2)
    filename = tempname()*".h5"
    BPJSpec.touch(filename,Nfreq,Nbase,Ntime)

    data  = rand(Complex64,Nbase,Nfreq)
    flags = zeros(Bool,Nbase,Nfreq)
    BPJSpec.grid(filename,0.0,data,flags)

    h5open(filename,"r") do file
        for β = 1:Nfreq
            group = file[string(β)]
            @test squeeze(group["data"][:,:,1],3) == reinterpret(Float32,data[:,β],(2,Nbase))
            @test squeeze(group["weights"][:,1],2) == ones(Nbase)
            @test group["data"][:,:,2:Ntime] == zeros(2,Nbase,Ntime-1)
            @test group["weights"][:,2:Ntime] == zeros(Nbase,Ntime-1)
        end
    end
end

