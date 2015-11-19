# Test NoiseModel
let
    Tsys0 = 6420
    α = 2.56
    ν0 = 45e6
    Δν = 24e3
    τ_total = 24*3600.0
    τ_int   = 30.0
    noise = NoiseModel(Tsys0,α,ν0,Δν,τ_total,τ_int)
    @test noise(0,ν0) ≈ Tsys0^2/(τ_total*Δν)
    @test noise(0,2ν0) ≈ Tsys0^2/(τ_total*Δν) * 2^(-2α)
end

# Test ForegroundModel
let
    ν0 = 45e6
    A = 123
    α = 1.0
    β = 1.0
    ζ = Inf
    foreground = ForegroundModel(ν0,A,α,β,ζ)
    @test foreground(0,ν0,ν0) == A
    @test foreground(20,ν0,ν0) / foreground(10,ν0,ν0) ≈ ((20+1)/(10+1))^(-α)
    @test foreground(0,2ν0,2ν0) / foreground(0,ν0,ν0) ≈ 2^(-2β)
end

# Test SignalModel
let
    kpara = linspace(0.5,1.0,1001)
    kperp = linspace(0.0,0.1,1000)
    power = ones(length(kpara),length(kperp))
    signal = SignalModel(kpara,kperp,power)

    ν = 45e6
    z = BPJSpec.redshift(ν)
    χ = BPJSpec.comoving_distance(z)

    ν1 = 45e6
    ν2 = 45e6 + 1e6
    z1 = BPJSpec.redshift(ν1)
    z2 = BPJSpec.redshift(ν2)
    χ1 = BPJSpec.comoving_distance(z1)
    χ2 = BPJSpec.comoving_distance(z2)
    Δχ = χ2 - χ1

    for l in (0, 50, 100)
        @test signal(l,ν,ν) ≈ (kpara[end]-kpara[1]) / (π*χ^2)
        @test signal(l,ν1,ν2) ≈ (sin(kpara[end]*Δχ) - sin(kpara[1]*Δχ)) / (π*χ1*χ2*Δχ)
        # test that the ν1 ≠ ν2 branch approaches the ν1 = ν2 branch
        # in the limit that ν2 → ν1
        @test signal(l,ν,ν) ≈ signal(l,ν,ν+1.0)
    end

    # ok, now let's try with a power spectrum that is evolving with k
    for j = 1:length(kperp), i = 1:length(kpara)
        power[i,j] = kpara[i] + kperp[j]
    end
    signal = SignalModel(kpara,kperp,power)

    for l in (0, 50, 100)
        P1 = kpara[  1] + l/χ
        P2 = kpara[end] + l/χ
        @test signal(l,ν,ν) ≈ 0.5 * (P1+P2) * (kpara[end]-kpara[1]) / (π*χ^2)
        P1 = kpara[  1] + 2l/(χ1+χ2)
        P2 = kpara[end] + 2l/(χ1+χ2)
        @test signal(l,ν1,ν2) ≈ (((P2*sin(kpara[end]*Δχ) - P1*sin(kpara[1]*Δχ)) / Δχ
                                 + (cos(kpara[end]*Δχ) - cos(kpara[1]*Δχ)) / Δχ^2) / (π*χ1*χ2))
    end
end

