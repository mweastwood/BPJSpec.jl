# Copyright (c) 2015-2017 Michael Eastwood
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

module FastTransformsWrapper

using FastTransforms
using CasaCore.Measures

struct SHT
    sph2fourier_plan
    #synthesis_plan
    #analysis_plan
    lmax :: Int
    mmax :: Int
    size :: Tuple{Int, Int}
end

function plan_sht(lmax, mmax, size)
    alm = zeros(lmax+1, 2mmax+1)
    map = zeros(size)
    sph2fourier_plan = plan_sph2fourier(alm, sketch=:none)
    #synthesis_plan = FastTransforms.plan_synthesis(map)
    #analysis_plan = FastTransforms.plan_analysis(map)
    SHT(sph2fourier_plan, lmax, mmax, size)
end

function plan_sht(metadata, size)
    plan_sht(sphericalharmonics.lmax, sphericalharmonics.mmax, size)
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

Base.getindex(map::Map, idx...) = map.matrix[idx...]
Base.setindex!(map::Map, value, idx...) = map.matrix[idx...] = value
Base.size(map::Map) = size(map.matrix)

function alm2map(sht, alm)
    # convert to bivariate Fourier series
    fourier = sht.sph2fourier_plan*alm.matrix

    # pad the Fourier series to the desired map size
    #padded_fourier = zeros(eltype(fourier), sht.size)
    #padded_fourier[1:sht.lmax+1, 1:2sht.mmax+1] = fourier
    padded_fourier = fourier

    # synthesize the map
    synthesis_plan = FastTransforms.plan_synthesis(padded_fourier)
    output = A_mul_B!(zero(padded_fourier), synthesis_plan, padded_fourier)
    Map(output)
end

function map2alm(sht, map)
    # analyze the map
    analysis_plan = FastTransforms.plan_analysis(map.matrix)
    fourier = A_mul_B!(zero(map.matrix), analysis_plan, map.matrix)

    # cut the Fourier series down to the right size (for the desired lmax, mmax)
    #cut_fourier = fourier[1:sht.lmax+1, 1:2sht.mmax+1]
    cut_fourier = fourier

    # convert to spherical harmonic coefficients
    output = sht.sph2fourier_plan\cut_fourier
    Alm(sht.lmax, sht.mmax, output)
end

Base.:*(sht::SHT, map::Map) = map2alm(sht, map)
Base.:\(sht::SHT, alm::Alm) = alm2map(sht, alm)

function index2vector(map::Map, idx, jdx)
    n, m = size(map)
    θ = π*(idx-0.5)/n
    ϕ = π - 2π*(jdx-1)/m
    s = sin(θ)
    x = s*cos(ϕ)
    y = s*sin(ϕ)
    z = cos(θ)
    Direction(dir"ITRF", x, y, z)
end

end

