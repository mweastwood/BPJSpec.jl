let
    # Tests against Ned Wright's cosmology calculator
    # (www.astro.ucla.edu/~wright/CosmoCalc.html)
    # - remember to set H0 = 69, OmegaM = 0.29, flat
    @test (BPJSpec.comoving_distance(0.1) - 424.8) < 0.1
    @test (BPJSpec.comoving_distance(10.) - 9689.5) < 0.1
    @test (BPJSpec.age(0.1) - 12.465) < 0.001
    @test (BPJSpec.age(10.) - 0.479) < 0.001
end

for i = 1:10
    z  = 100*rand()
    ν  = BPJSpec.frequency(z)
    z′ = BPJSpec.redshift(ν)
    @test z ≈ z′
end

let
    @test BPJSpec.beam_solid_angle(SineBeam(0)) ≈ 2π
end

