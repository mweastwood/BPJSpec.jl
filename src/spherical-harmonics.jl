using FastTransforms
import LibHealpix

struct SHT
    sph2fourier_plan
    synthesis_plan
    analysis_plan
    lmax :: Int
    mmax :: Int
    size :: Tuple{Int, Int}
end

function plan_sht(lmax, mmax, size)
    alm = zeros(lmax+1, 2mmax+1)
    map = zeros(size)
    sph2fourier_plan = plan_sph2fourier(alm)
    synthesis_plan = FastTransforms.plan_synthesis(map)
    analysis_plan = FastTransforms.plan_analysis(map)
    SHT(sph2fourier_plan, synthesis_plan, analysis_plan, lmax, mmax, size)
end

function map2alm(map, lmax, mmax)
    # analyze the map
    analysis_plan = FastTransforms.plan_analysis(map.matrix)
    fourier = A_mul_B!(zero(map.matrix), analysis_plan, map.matrix)

    # convert to spherical harmonic coefficients
    plan = plan_sph2fourier(fourier)
    output = plan\fourier
    Alm(lmax, mmax, output)
end

struct Alm
    lmax :: Int
    mmax :: Int
    matrix :: Matrix{Float64}
end

function Alm(lmax, mmax)
    matrix = zeros(lmax+1, 2mmax+1)
    Alm(lmax, mmax, matrix)
end

function Base.getindex(alm::Alm, l, m)
    idx = l - m + 1
    jdx = 2m + 1
    if m == 0
        return complex(alm.matrix[idx, jdx])
    else
        return complex(alm.matrix[idx, jdx], alm.matrix[idx, jdx-1]) / √2
    end
end

function Base.setindex!(alm::Alm, value, l, m)
    idx = l - m + 1
    jdx = 2m + 1
    if m == 0
        alm.matrix[idx, jdx] = real(value)
        return complex(real(value))
    else
        sqrt2 = √2
        alm.matrix[idx, jdx]   = real(value) * sqrt2
        alm.matrix[idx, jdx-1] = imag(value) * sqrt2
        return value
    end
end

struct Map <: AbstractMatrix{Float64}
    matrix :: Matrix{Float64}
end

Base.getindex(map::Map, idx) = map.matrix[idx]
Base.size(map::Map) = size(map.matrix)

function alm2map(sht, alm)
    # convert to bivariate Fourier series
    @time fourier = sht.sph2fourier_plan*alm.matrix

    # pad the Fourier series to the desired map size
    padded_fourier = zeros(eltype(fourier), sht.size)
    padded_fourier[1:sht.lmax+1, 1:2sht.mmax+1] = fourier

    # synthesize the map
    @time output = A_mul_B!(zero(padded_fourier), sht.synthesis_plan, padded_fourier)
    Map(output)
end

function map2alm(sht, map)
    # analyze the map
    @time fourier = A_mul_B!(zero(map.matrix), sht.analysis_plan, map.matrix)

    # cut the Fourier series down to the right size (for the desired lmax, mmax)
    cut_fourier = fourier[1:sht.lmax+1, 1:2sht.mmax+1]

    # convert to spherical harmonic coefficients
    @time output = sht.sph2fourier_plan\cut_fourier
    Alm(sht.lmax, sht.mmax, output)
end

Base.:*(sht::SHT, map::Map) = map2alm(sht, map)
Base.:\(sht::SHT, alm::Alm) = alm2map(sht, alm)

#using PyPlot

function test(lmax, mmax)
    alm1 = LibHealpix.Alm(Complex128, lmax, mmax)
    for m = 0:mmax, l = m:lmax
        LibHealpix.@lm alm1[l, m] = complex(randn(), randn())
    end

    alm2 = Alm(lmax, mmax)
    for m = 0:mmax, l = m:lmax
        alm2[l, m] = LibHealpix.@lm alm1[l, m]
    end

    @time sht = plan_sht(lmax, mmax, (2048, 4095))

    @time healpix = LibHealpix.alm2map(alm1, 2048)
    println("---")
    @time map2 = sht\alm2
    println("---")

    @show vecnorm(map2-map1)/vecnorm(map1)

    map1 = similar(map2)
    n, m = size(map2)
    θ = (0.5:n-0.5)*π/n
    ϕ = (0:m-1)*2π/m
    for jdx = 1:m, idx = 1:n
        map1[idx, jdx] = LibHealpix.interpolate(healpix, θ[idx], π-ϕ[jdx])
    end

    @time alm1′ = LibHealpix.map2alm(healpix, lmax, mmax, iterations=2)
    println("---")
    @time alm2′ = sht*map2
    println("---")
    @show vecnorm(alm1-alm1′)/vecnorm(alm1)
    @show vecnorm(alm2.matrix-alm2′.matrix)/vecnorm(alm2.matrix)


    #figure(1); clf()
    #subplot(3, 1, 1)
    #imshow(map1)
    #colorbar()
    #subplot(3, 1, 2)
    #imshow(map2)
    #colorbar()
    #subplot(3, 1, 3)
    #imshow(log10.(abs.(map2.-map1)))
    #colorbar()
    #tight_layout()

    nothing
end

