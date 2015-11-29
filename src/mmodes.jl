# Copyright (c) 2015 Michael Eastwood
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

immutable MModesMeta <: Metadata
    m::UnitRange{Int}
    ν::Vector{Float64}

    function MModesMeta(m,ν)
        if length(m) > 1 && length(ν) > 1
            error("Cannot simultaneously have multiple values of m and multiple frequency channels.")
        end
        new(m,ν)
    end
end

MModesMeta(m::Int,ν::AbstractVector) = MModesMeta(m:m,collect(ν))
MModesMeta(mmax::Int,ν::Float64) = MModesMeta(0:mmax,[ν])

==(lhs::MModesMeta,rhs::MModesMeta) = lhs.m == rhs.m && lhs.ν == rhs.ν

doc"""
    typealias MModes Blocks{VectorBlock, MModesMeta}

This type stores the $m$-modes organized into blocks.

$m$-modes are the Fourier transform of a visibility with
respect to time. These are relate to the spherical harmonic
coefficients of the sky brightness through a matrix equation.

\\[
    v = Ba,
\\]

where $v$ is the vector of $m$-modes, $B$ is the transfer
matrix, and $a$ is the vector of spherical harmonic coefficients.
"""
typealias MModes Blocks{VectorBlock, MModesMeta}

initial_block_size(::Type{MModes}, Nbase, m) = (two(m)*Nbase,)

function call(::Type{MModes}, Nbase::Int, mmax::Int, ν::Float64)
    meta = MModesMeta(mmax,ν)
    blocks = VectorBlock[]
    for m = 0:mmax
        sz = initial_block_size(MModes,Nbase,m)
        push!(blocks,VectorBlock(sz))
    end
    MModes(blocks,meta)
end

function call(::Type{MModes}, vector::Vector{Complex128}, m::Int, ν::AbstractVector)
    # When foreground-filtering the m-modes, it is often convenient to call
    # `full(::MModes)` to turn the m-modes into a standard vector. This function
    # is the inverse of `full`.
    Nfreq = length(ν)
    N = div(length(vector),Nfreq)
    blocks = [VectorBlock(vector[(β-1)*N+1:β*N]) for β = 1:Nfreq]
    meta = MModesMeta(m:m,collect(ν))
    MModes(blocks,meta)
end

is_single_frequency(meta::MModesMeta) = length(meta.ν) == 1
is_single_m(meta::MModesMeta) = length(meta.m) == 1
is_single_frequency(v::MModes) = is_single_frequency(v.meta)
is_single_m(v::MModes) = is_single_m(v.meta)

mmax(meta::MModesMeta) = maximum(meta.m)
mmax(v::MModes) = mmax(v.meta)

Nfreq(meta::MModesMeta) = length(meta.ν)
Nfreq(v::MModes) = Nfreq(v.meta)

"""
    mmodes(visibilities; frequency=0.0, mmax=100)

Calculate the m-modes from the given visibilities. The visibilities
should be provided as a matrix where the first dimension indicates
the baseline and the second dimension indicates the sidereal time.
The FFT will be performed over the second dimension. The
visibilities should span a full sidereal day.
"""
function mmodes{T<:Complex}(visibilities::Matrix{T};
                            mmax::Int = 100,
                            frequency::Float64 = 0.0)
    Nbase,Ntime = size(visibilities)
    M = fft(visibilities,2)/Ntime
    v = MModes(Nbase,mmax,frequency)
    for α = 1:Nbase
        v[1][α] = M[α,1]
    end
    for m = 1:mmax, α = 1:Nbase
        α1 = α         # positive m
        α2 = α + Nbase # negative m
        v[m+1][α1] =      M[α,m+1]
        v[m+1][α2] = conj(M[α,Ntime+1-m])
    end
    v
end

"""
    visibilities(v::MModes)

Calculate the visibilities from the given m-modes. The visibilities
will be returned as a two dimensional array where the first dimension
indicates the baseline and the second dimension indicates the sidereal
time, which evenly spans a full sidereal day.
"""
function visibilities(v::MModes)
    is_single_frequency(v) || error("Expected single-frequency m-modes.")
    Nbase = size(v[1].block,1)
    Ntime = 2mmax(v)+1
    M = zeros(Complex128,Nbase,Ntime)
    for α = 1:Nbase
        M[α,1] = v[1][α]
    end
    for m = 1:mmax(v), α = 1:Nbase
        α1 = α         # positive m
        α2 = α + Nbase # negative m
        M[α,m+1]       =      v[m+1][α1]
        M[α,Ntime+1-m] = conj(v[m+1][α2])
    end
    ifft(M,2)*Ntime
end

function save(filename, v::MModes)
    if !isfile(filename)
        jldopen(filename,"w",compress=true) do file
            file["description"] = "m-modes"
        end
    end

    jldopen(filename,"r+",compress=true) do file
        if read(file["description"]) != "m-modes"
            error("Attempting to write to a file that does not contain m-modes.")
        end

        if is_single_frequency(v)
            ν = v.meta.ν[1]
            name = @sprintf("%.3fMHz",ν/1e6)
            name in names(file) || g_create(file,name)
            group = file[name]
            for m = 0:mmax(v)
                block = v[m+1]
                group[string(m)] = block.block
            end
        elseif is_single_m(v)
            m = v.meta.m[1]
            for β = 1:Nfreq(v)
                ν = v.meta.ν[β]
                name = @sprintf("%.3fMHz",ν/1e6)
                name in names(file) || g_create(file,name)
                group = file[name]
                block = v[β]
                group[string(m)] = block.block
            end
        end
    end
end

function load(filename, meta::MModesMeta)
    blocks = VectorBlock[]
    jldopen(filename,"r") do file
        if read(file["description"]) != "m-modes"
            error("Attempting to read from a file that does not contain m-modes.")
        end

        if is_single_frequency(meta)
            ν = meta.ν[1]
            name = @sprintf("%.3fMHz",ν/1e6)
            group = file[name]
            for m = 0:mmax(meta)
                block = group[string(m)] |> read
                push!(blocks,VectorBlock(block))
            end
        elseif is_single_m(meta)
            m = meta.m[1]
            for β = 1:Nfreq(meta)
                ν = meta.ν[β]
                name = @sprintf("%.3fMHz",ν/1e6)
                group = file[name]
                block = group[string(m)] |> read
                push!(blocks,VectorBlock(block))
            end
        end
    end
    MModes(blocks,meta)
end

