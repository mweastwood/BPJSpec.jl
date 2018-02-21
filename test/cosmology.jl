@testset "cosmology.jl" begin
    # Tests against Ned Wright's cosmology calculator
    # (www.astro.ucla.edu/~wright/CosmoCalc.html)
    # - remember to set H0 = 69, OmegaM = 0.29, flat
    @test abs(BPJSpec.comoving_distance(0.1) -  424.8*u"Mpc") < 0.1*u"Mpc"
    @test abs(BPJSpec.comoving_distance(10.) - 9689.5*u"Mpc") < 0.1*u"Mpc"
    @test abs(BPJSpec.age(0.1) - 12.465*u"Gyr") < 0.001*u"Gyr"
    @test abs(BPJSpec.age(10.) -  0.479*u"Gyr") < 0.001*u"Gyr"

    for z in collect(0:1.23:30)
        ν  = BPJSpec.frequency(z)
        z′ = BPJSpec.redshift(ν)
        @test z ≈ z′
    end
end

