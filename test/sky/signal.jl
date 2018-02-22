@testset "signal.jl" begin
    kpara = collect(linspace(0.0, 1.0, 1001)) .* u"Mpc^-1"
    kperp = collect(linspace(0.0, 0.1, 1000)) .* u"Mpc^-1"
    power = ones(length(kpara), length(kperp)) .* u"mK^2*Mpc^3"
    signal = BPJSpec.SignalModel((10, 30), kpara, kperp, power)
    units = u"mK^2*Mpc^3"

    ν = 45u"MHz"
    z = BPJSpec.redshift(ν)
    χ = BPJSpec.comoving_distance(z)

    ν1 = 45u"MHz"
    ν2 = ν1 + 1u"MHz"
    z1 = BPJSpec.redshift(ν1)
    z2 = BPJSpec.redshift(ν2)
    χ1 = BPJSpec.comoving_distance(z1)
    χ2 = BPJSpec.comoving_distance(z2)
    Δχ = χ2 - χ1

    for l in (0, 50, 100)
        @test signal(l, ν, ν) ≈ (kpara[end]-kpara[1]) / (π*χ^2) * units
        @test signal(l, ν1, ν2) ≈ (sin(kpara[end]*Δχ) - sin(kpara[1]*Δχ)) / (π*χ1*χ2*Δχ) * units
        # test that the ν1 ≠ ν2 branch approaches the ν1 = ν2 branch
        # in the limit that ν2 → ν1
        @test signal(l, ν, ν) ≈ signal(l, ν, ν+1u"Hz")
    end

    # ok, now let's try with a power spectrum that is evolving with k
    for j = 1:length(kperp), i = 1:length(kpara)
        power[i, j] = ustrip(kpara[i] + kperp[j])*units
    end
    signal = BPJSpec.SignalModel((10, 30), kpara, kperp, power)

    for l in (0, 50, 100)
        P1 = ustrip(kpara[  1] + l/χ)*units
        P2 = ustrip(kpara[end] + l/χ)*units
        @test signal(l, ν, ν) ≈ 0.5 * (P1+P2) * (kpara[end]-kpara[1]) / (π*χ^2)
        P1 = ustrip(kpara[  1] + 2l/(χ1+χ2))*units
        P2 = ustrip(kpara[end] + 2l/(χ1+χ2))*units
        @test signal(l, ν1, ν2) ≈ (((P2*sin(kpara[end]*Δχ) - P1*sin(kpara[1]*Δχ)) / Δχ
                                    + (cos(kpara[end]*Δχ) - cos(kpara[1]*Δχ)) / Δχ^2 * units*u"Mpc")
                                    / (π*χ1*χ2))
    end
end
