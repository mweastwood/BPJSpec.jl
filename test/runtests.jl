using BPJSpec
if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end
using CasaCore.Measures
using CasaCore.Tables
using LibHealpix
using TTCal
using JLD

include("setup.jl")

srand(123)
@testset "BPJSpec Tests" begin
    include("special.jl")
    include("physics.jl")
    include("visibilities.jl")
    #include("mmodes.jl")
    #include("transfermatrix.jl")
    #include("alm.jl")

    #include("noise.jl")
    #include("sky.jl")
end

