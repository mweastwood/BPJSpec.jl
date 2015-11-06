using BPJSpec
using Base.Test
using CasaCore.Measures
using CasaCore.Tables
using TTCal
using JLD

function xyz2uvw(x,y,z)
    Nant = length(x)
    Nbase = div(Nant*(Nant-1),2) + Nant
    u = Array{Float64}(Nbase)
    v = Array{Float64}(Nbase)
    w = Array{Float64}(Nbase)
    α = 1
    for i = 1:Nant, j = i:Nant
        u[α] = x[j]-x[i]
        v[α] = y[j]-y[i]
        w[α] = z[j]-z[i]
        α += 1
    end
    u,v,w
end

function ant1ant2(Nant)
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1 = Array(Int32,Nbase)
    ant2 = Array(Int32,Nbase)
    α = 1
    for i = 1:Nant, j = i:Nant
        ant1[α] = i
        ant2[α] = j
        α += 1
    end
    ant1,ant2
end

function createms(Nant,Nfreq)
    Nbase = div(Nant*(Nant-1),2) + Nant

    x = 100*randn(Nant)
    y = 100*randn(Nant)
    z = randn(Nant)
    u,v,w = xyz2uvw(x,y,z)
    ν = linspace(40e6,60e6,Nfreq) |> collect
    t = (2015.-1858.)*365.*24.*60.*60. # a rough current Julian date (in seconds)
    ant1,ant2 = ant1ant2(Nant)

    frame = ReferenceFrame()
    pos = observatory("OVRO_MMA")
    set!(frame,Epoch(epoch"UTC",Quantity(t,"s")))
    set!(frame,pos)
    zenith = Direction(dir"AZEL",q"0.0deg",q"90.0deg")
    phase_dir = measure(frame,zenith,dir"J2000")

    name  = tempname()*".ms"
    table = Table(name)

    subtable = Table("$name/SPECTRAL_WINDOW")
    Tables.addRows!(subtable,1)
    subtable["CHAN_FREQ"] = reshape(ν,length(ν),1)
    unlock(subtable)

    subtable = Table("$name/ANTENNA")
    Tables.addRows!(subtable,Nant)
    x,y,z = Measures.xyz_in_meters(pos)
    subtable["POSITION"] = [x;y;z]*ones(1,Nant)
    unlock(subtable)

    subtable = Table("$name/FIELD")
    Tables.addRows!(subtable,1)
    subtable["PHASE_DIR"] = reshape([longitude(phase_dir);latitude(phase_dir)],2,1)
    unlock(subtable)

    Tables.addRows!(table,Nbase)
    table[kw"SPECTRAL_WINDOW"] = "Table: $name/SPECTRAL_WINDOW"
    table[kw"ANTENNA"] = "Table: $name/ANTENNA"
    table[kw"FIELD"] = "Table: $name/FIELD"
    table["ANTENNA1"] = ant1-1
    table["ANTENNA2"] = ant2-1
    table["UVW"] = [u v w]'
    table["TIME"] = fill(float(t),Nbase)
    table["FLAG_ROW"] = zeros(Bool,Nbase)
    table["FLAG"] = zeros(Bool,4,Nfreq,Nbase)
    table["DATA"] = zeros(Complex64,4,Nfreq,Nbase)
    unlock(table)

    name,MeasurementSet(name)
end

srand(123)
include("cornerturn.jl")

#=
# Test special functions
for i = 1:10
    θ =  π*rand()
    ϕ = 2π*rand()
    @test BPJSpec.Y(0,0,θ,ϕ) ≈ 1/2 * sqrt(1/π)
    @test BPJSpec.Y(1,0,θ,ϕ) ≈ 1/2 * sqrt(3/π) * cos(θ)
    @test BPJSpec.Y(1,+1,θ,ϕ) ≈ -1/2 * sqrt(3/(2π)) * sin(θ) * exp(+1im*ϕ)
    @test BPJSpec.Y(1,-1,θ,ϕ) ≈ +1/2 * sqrt(3/(2π)) * sin(θ) * exp(-1im*ϕ)
end

# Test the plane wave expansion
using LibHealpix
let u = 2.0, v = -3.0, w = 2.0
    alm_real, alm_imag = BPJSpec.planewave(u,v,w)
    map_real = alm2map(alm_real,nside=512)
    map_imag = alm2map(alm_imag,nside=512)
    map_test_real = HealpixMap(Float64,512)
    map_test_imag = HealpixMap(Float64,512)
    for i = 1:nside2npix(512)
        vec = LibHealpix.pix2vec_ring(512,i)
        phase = 2π*(u*vec[1]+v*vec[2]+w*vec[3])
        map_test_real[i] = cos(phase)
        map_test_imag[i] = sin(phase)
    end
    @test pixels(map_real) ≈ pixels(map_test_real)
    @test pixels(map_imag) ≈ pixels(map_test_imag)
end

# Test the matrix trace
for i = 1:10
    A = rand(Complex128,50,50)
    B = rand(Complex128,50,50)
    @test trace(A*B) ≈ BPJSpec.tr(A,B)
end

# Test force_hermitian and force_posdef
let
    A = rand(10,10)
    @test ishermitian(BPJSpec.force_hermitian(A))
    B = diagm(linspace(1,-1e-11,10))
    @test isposdef(BPJSpec.force_posdef(B))
end

# Tests against Ned Wright's cosmology calculator
# (www.astro.ucla.edu/~wright/CosmoCalc.html)
# - remember to set H0 = 69, OmegaM = 0.29, flat
@test (BPJSpec.comoving_distance(0.1) - 424.8) < 0.1
@test (BPJSpec.comoving_distance(10.) - 9689.5) < 0.1
@test (BPJSpec.age(0.1) - 12.465) < 0.001
@test (BPJSpec.age(10.) - 0.479) < 0.001

# Test visibility generation
let Nant = 3, mmax = 5
    Nbase = div(Nant*(Nant-1),2)
    data = zeros(Complex128,Nbase,2mmax+1)
    rand!(data)
    mmodes = MModes(data,mmax=mmax)
    data′ = visibilities(mmodes)
    @test data ≈ data′
end

# Test ForegroundModel
let
    ν0 = 45.0
    A = 123
    α = 1.0
    β = 1.0
    ζ = 1.0
    foreground = ForegroundModel(ν0,A,α,β,ζ)
    @test foreground(0,ν0,ν0) == A
end

# Test SphericalSignalModel
let
    k = logspace(log10(0.03),log10(0.3))
    P = ones(length(k)-1)
    signal = SphericalSignalModel(k,P)

    ν = 45.0
    z = BPJSpec.redshift(ν)
    r = BPJSpec.comoving_distance(z)
    @test BPJSpec.Csignal_spherical(0,ν,ν,k,P) ≈ (k[end]-k[1])/(π*r^2)
    @test signal(0,ν,ν) ≈ (k[end]-k[1])/(π*r^2)

    ν1 = 45.0
    ν2 = 45.1
    z1 = BPJSpec.redshift(ν1)
    z2 = BPJSpec.redshift(ν2)
    r1 = BPJSpec.comoving_distance(z1)
    r2 = BPJSpec.comoving_distance(z2)
    Δr = r2-r1
    @test(BPJSpec.Csignal_spherical(0,ν1,ν2,k,P)
            ≈ (sin(k[end]*Δr)-sin(k[1]*Δr))/(π*Δr*r1*r2))
    @test signal(0,ν1,ν2) ≈ (sin(k[end]*Δr)-sin(k[1]*Δr))/(π*Δr*r1*r2)
end

# Test SpectralMModes
let
    mmax = 10
    a = BPJSpec.MModeBlock{mmax}(0,complex(randn(5),randn(5)))
    b = BPJSpec.MModeBlock{mmax}(0,complex(randn(5),randn(5)))
    v = SpectralMModes{mmax}(0,[a,b])
    @test full(v) == [a.block; b.block]
end

# Test SpectralTransferMatrix
let
    lmax = 10
    mmax = 10
    a = BPJSpec.TransferMatrixBlock{lmax,mmax}(0,complex(randn(5,5),randn(5,5)))
    b = BPJSpec.TransferMatrixBlock{lmax,mmax}(0,complex(randn(5,5),randn(5,5)))
    z = zeros(Complex128,5,5)
    B = SpectralTransferMatrix{lmax,mmax}(0,[a,b])
    @test full(B) == [a.block z;
                      z b.block]
end
=#

