using FastTransforms
import LibHealpix

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

function alm2map(alm)
    # convert to bivariate Fourier series
    @time plan = SlowSphericalHarmonicPlan(alm.matrix)
    @time fourier = plan*alm.matrix

    # pad the Fourier series to the desired map size
    #N, M = size(fourier)
    #padded_fourier = zeros(eltype(fourier), sz)
    #@views padded_fourier[1:N, 1:M] = fourier

    # synthesize the map
    @time synthesis_plan = FastTransforms.plan_synthesis(fourier)
    @show typeof(synthesis_plan)
    @time output = A_mul_B!(zero(fourier), synthesis_plan, fourier)
    Map(output)
end

function map2alm(alm)
    # convert to bivariate Fourier series
    @time plan = SlowSphericalHarmonicPlan(alm.matrix)
    @time fourier = plan*alm.matrix

    # pad the Fourier series to the desired map size
    #N, M = size(fourier)
    #padded_fourier = zeros(eltype(fourier), sz)
    #@views padded_fourier[1:N, 1:M] = fourier

    # synthesize the map
    @time synthesis_plan = FastTransforms.plan_synthesis(fourier)
    @time output = A_mul_B!(zero(fourier), synthesis_plan, fourier)
    Map(output)
end


using PyPlot

function test(lmax, mmax)
    alm1 = LibHealpix.Alm(Complex128, lmax, mmax)
    for m = 0:mmax, l = m:lmax
        LibHealpix.@lm alm1[l, m] = complex(randn(), randn())
    end

    alm2 = Alm(lmax, mmax)
    for m = 0:mmax, l = m:lmax
        alm2[l, m] = LibHealpix.@lm alm1[l, m]
    end

    @time healpix = LibHealpix.alm2map(alm1, 2048)
    @time map2 = alm2map(alm2)

    map1 = similar(map2)
    n, m = size(map2)
    θ = (0.5:n-0.5)*π/n
    ϕ = (0:m-1)*2π/m
    for jdx = 1:m, idx = 1:n
        map1[idx, jdx] = LibHealpix.interpolate(healpix, θ[idx], π-ϕ[jdx])
    end

    figure(1); clf()
    subplot(3, 1, 1)
    imshow(map1)
    colorbar()
    subplot(3, 1, 2)
    imshow(map2)
    colorbar()
    subplot(3, 1, 3)
    imshow(log10.(abs.(map2.-map1)))
    colorbar()
    tight_layout()

    vecnorm(map2-map1)/vecnorm(map1)
end

