@testset "FastTransformsWrapper.jl" begin
    lmax = mmax = 5
    size = (lmax+1, 2mmax+1)
    plan = BPJSpec.FastTransformsWrapper.plan_sht(lmax, mmax, size)

    alm = BPJSpec.FastTransformsWrapper.Alm(lmax, mmax)
    alm[0, 0] = 1
    map = plan \ alm
    @test all(map .≈ 1/sqrt(4π))
    alm′ = plan * map
    @test alm′[0, 0] ≈ 1
end

