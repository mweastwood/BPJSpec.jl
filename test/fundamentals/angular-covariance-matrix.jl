# NOTE: It's very important that this model is linear in frequency because the tests rely on the
# fact that the bandwidth smeared version of this is the same as without bandwidth smearing.
struct SillySkyComponent <: BPJSpec.SkyComponent end
(::SillySkyComponent)(l, ν1, ν2) = (l + ustrip(uconvert(u"MHz", ν1 + ν2))) * u"K^2"

# NOTE: This model tries to engineer a situation where the bandwidth smearing forces the signal to
# zero.
struct CoolSkyComponent <: BPJSpec.SkyComponent end
function (::CoolSkyComponent)(l, ν1, ν2)
    x = ustrip(uconvert(u"MHz", ν1))/100 - 1
    if x ≥ 0
        return (1-20x) * u"K^2"
    else
        return (1+20x) * u"K^2"
    end
end

@testset "angular-covariance-matrix.jl" begin
    lmax = 2
    frequencies = [45u"MHz", 74u"MHz"]
    bandwidth   = [ 1u"Hz",   2u"Hz" ]
    component = BPJSpec.galactic_synchrotron()
    C = create(AngularCovarianceMatrix, NoFile(), component, lmax, frequencies, bandwidth)
    @test C[L(1), 0] == C[L(1), 1] # independent of m
    for l = 0:lmax
        Cl = C[L(l), 0]
        @test ishermitian(Cl)
        @test isposdef(Cl)
        for β1 = 1:length(frequencies), β2 = 1:length(frequencies)
            ν1 = frequencies[β1]
            ν2 = frequencies[β2]
            # approximate because we haven't accounted for the bandwidth smearing (just made it
            # small by setting the bandwidth to be small)
            @test Cl[β1, β2]*u"K^2" ≈ component(l, ν1, ν2)
        end
    end

    # Now we'd really like to test that the bandwidth smearing routine is working correctly, so
    # let's verify that with simpler models that we can check analytically.

    component = SillySkyComponent()
    frequencies = [45u"MHz", 74u"MHz"]
    bandwidth   = [20u"MHz", 20u"MHz"] # add a lot of bandwidth to make bandwidth smearing important
    C = create(AngularCovarianceMatrix, NoFile(), component, lmax, frequencies, bandwidth)
    for l = 0:lmax
        Cl = C[L(l), 0]
        for β1 = 1:length(frequencies), β2 = 1:length(frequencies)
            ν1 = frequencies[β1]
            ν2 = frequencies[β2]
            Δν1 = bandwidth[β1]
            Δν2 = bandwidth[β2]
            @test Cl[β1, β2]*u"K^2" ≈ component(l, ν1, ν2)
        end
    end

    component = CoolSkyComponent()
    # carefully chosen frequency and bandwidth to make this smear to zero
    frequencies = [100u"MHz"]
    bandwidth   = [20u"MHz"]
    C = create(AngularCovarianceMatrix, NoFile(), component, lmax, frequencies, bandwidth)
    for l = 0:lmax
        Cl = C[L(l), 0]
        for β1 = 1:length(frequencies), β2 = 1:length(frequencies)
            ν1 = frequencies[β1]
            ν2 = frequencies[β2]
            Δν1 = bandwidth[β1]
            Δν2 = bandwidth[β2]
            # note that the tolerance here is controlled by the tolerance used for cubature
            @test abs(Cl[β1, β2]*u"K^2") < 1e-8u"K^2"
        end
    end
end

