using BPJSpec
using Base.Test

srand(123)

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

let Nant = 3, mmax = 5
    Nbase = div(Nant*(Nant-1),2)
    data = zeros(Complex128,Nbase,2mmax+1)
    rand!(data)
    mmodes = MModes(data,mmax=mmax)
    data′ = visibilities(mmodes)
    @test data ≈ data′
end

