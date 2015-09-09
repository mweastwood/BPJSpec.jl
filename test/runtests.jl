using BPJSpec
using Base.Test

srand(123)

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
using HEALPix
let u = 2.0, v = -3.0, w = 2.0
    alm_real, alm_imag = BPJSpec.planewave(u,v,w)
    map_real = alm2map(alm_real,nside=512)
    map_imag = alm2map(alm_imag,nside=512)
    map_test_real = HEALPixMap(Float64,512)
    map_test_imag = HEALPixMap(Float64,512)
    for i = 1:nside2npix(512)
        vec = HEALPix.pix2vec_ring(512,i)
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

let
    ν = 45.0
    z = BPJSpec.redshift(ν)
    r = BPJSpec.comoving_distance(z)
    k = logspace(log10(0.03),log10(0.3))
    P = ones(length(k)-1)
    @test BPJSpec.Csignal_spherical(0,ν,ν,k,P) ≈ (k[end]-k[1])/(π*r^2)
end

