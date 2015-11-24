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

doc"""
    immutable MModeBlock <: AbstractVectorBlock

This type stores the $m$-modes corresponding to a single block of
the transfer matrix.

**Fields:**

* `blocks` stores the actual $m$-modes
* `m` stores the value of $m$ (the azimuthal spherical harmonic quantum number)
* `ν` stores the frequency (in Hz)
"""
immutable MModeBlock <: AbstractVectorBlock
    block::Vector{Complex128}
    m::Int
    ν::Float64
end

default_size(::Type{MModeBlock},Nbase,m) = (two(m)*Nbase,)
metadata(v::MModeBlock) = (v.m,v.ν)

doc"""
    MModeBlock(Nbase,m,ν)

Create an empty block of $m$-modes of the default size.
"""
function MModeBlock(Nbase,m,ν)
    sz = default_size(MModeBlock,Nbase,m)
    block = zeros(Complex128,sz)
    MModeBlock(block,m,ν)
end

doc"""
    immutable MModes{rep} <: AbstractBlockVector

This type represents a collection of $m$-modes.

$m$-modes are the Fourier transform of a visibility with
respect to time. These are relate to the spherical harmonic
coefficients of the sky brightness through a matrix equation.

\\[
    v = Ba,
\\]

where $v$ is the vector of $m$-modes, $B$ is the transfer
matrix, and $a$ is the vector of spherical harmonic coefficients.

See the documentation for `TransferMatrix` for a description
of the type parameter.

**Fields:**

* `blocks` is a list of `MModeBlock`s
"""
immutable MModes{rep} <: AbstractBlockVector
    blocks::Vector{MModeBlock}
end

doc"""
    MModes(blocks)

Construct a vector of $m$-modes from the list of blocks.

The representation (`one_ν` or `one_m`) is inferred from the provided
list of blocks.
"""
function MModes(blocks)
    m = [block.m for block in blocks]
    ν = [block.ν for block in blocks]
    mmax = length(m) - 1

    # if all the blocks are the same frequency and are sorted in order
    # of increasing m, we have the one_ν representation
    if length(unique(ν)) == 1 && m == collect(0:mmax)
        return MModes{one_ν}(blocks)
    end

    # if all the blocks have different frequencies (not necessarily sorted)
    # but only have one value of m, we have the one_m representation
    if length(unique(ν)) == length(ν) && length(unique(m)) == 1
        return MModes{one_m}(blocks)
    end

    error("""
    The list of blocks must either:

    * have a consistent frequency channel representing all values of m from 0 to some mmax, or
    * have different frequency channels, but a consistent value of m
    """)
end

function MModes(Nbase::Int,mmax::Int,ν::Float64)
    [MModeBlock(Nbase,m,ν) for m = 0:mmax] |> MModes
end

function MModes(vector::Vector{Complex128},m::Int,ν::AbstractVector)
    Nfreq = length(ν)
    N = div(length(vector),Nfreq)
    blocks = [MModeBlock(vector[(β-1)*N+1:β*N],m,ν[β]) for β = 1:Nfreq]
    MModes(blocks)
end

# α labels the baseline
# β labels the frequency channel
# m labels the mode

getindex(v::MModeBlock,α) = v.block[α]
setindex!(v::MModeBlock,x,α) = v.block[α] = x

getindex(v::MModes{one_ν},m) = v.blocks[m+1]
getindex(v::MModes{one_ν},α,m) = v[m][α]
setindex!(v::MModes{one_ν},x,α,m) = v[m][α] = x

getindex(v::MModes{one_m},β) = v.blocks[β]
getindex(v::MModes{one_m},α,β) = v[β][α]
setindex!(v::MModes{one_m},x,α,β) = v[β][α] = x

frequency(v::MModes{one_ν}) = v[0].ν
mmax(v::MModes{one_ν}) = length(v.blocks)-1

getm(v::MModes{one_m}) = v[1].m
frequencies(v::MModes{one_m}) = [v[β].ν for β = 1:length(v.blocks)]

"""
    MModes(visibilities, ν; mmax=100)

Calculate the m-modes from the given visibilities. The visibilities
should be provided as a matrix where the first dimension indicates
the baseline and the second dimension indicates the sidereal time.
The FFT will be performed over the second dimension. The
visibilities should span a full sidereal day.
"""
function MModes{T<:Complex}(visibilities::Matrix{T},ν;
                            mmax::Int = 100)
    Nbase,Ntime = size(visibilities)
    M = fft(visibilities,2)/Ntime
    v = MModes(Nbase,mmax,ν)
    for α = 1:Nbase
        v[α,0] = M[α,1]
    end
    for m = 1:mmax, α = 1:Nbase
        α1 = α         # positive m
        α2 = α + Nbase # negative m
        v[α1,m] =      M[α,m+1]
        v[α2,m] = conj(M[α,Ntime+1-m])
    end
    v
end

"""
    visibilities(v::MModes{one_ν})

Calculate the visibilities from the given m-modes. The visibilities
will be returned as a two dimensional array where the first dimension
indicates the baseline and the second dimension indicates the sidereal
time, which evenly spans a full sidereal day.
"""
function visibilities(v::MModes{one_ν})
    Nbase = size(v[0].block,1)
    Ntime = 2mmax(v)+1
    M = zeros(Complex128,Nbase,Ntime)
    for α = 1:Nbase
        M[α,1] = v[α,0]
    end
    for m = 1:mmax(v), α = 1:Nbase
        α1 = α         # positive m
        α2 = α + Nbase # negative m
        M[α,m+1]       =      v[α1,m]
        M[α,Ntime+1-m] = conj(v[α2,m])
    end
    ifft(M,2)*Ntime
end

function save_mmodes(filename, v::MModes{one_ν})
    if !isfile(filename)
        jldopen(filename,"w",compress=true) do file
            file["mmax"] = mmax(v)
        end
    end

    jldopen(filename,"r+",compress=true) do file
        read(file["mmax"]) == mmax(v) || error("mmax is inconsistent")

        name_ν  = @sprintf("%.3fMHz",frequency(v)/1e6)
        name_ν in names(file) && error("group $(name_ν) already exists in file")
        group_ν = g_create(file,name_ν)

        for m = 0:mmax(v)
            group_m = g_create(group_ν,string(m))
            group_m["block"] = v[m].block
        end
    end
end

function save_mmodes(filename, v::MModes{one_m})
    if !isfile(filename)
        jldopen(filename,"w",compress=true) do file
            file["mmax"] = getm(v)
        end
    end

    jldopen(filename,"r+",compress=true) do file
        m = getm(v)
        mmax = read(file["mmax"])
        if m > mmax
            if "mmax" in names(file)
                o_delete(file,"mmax")
            end
            file["mmax"] = m
        end

        ν = frequencies(v)
        for β = 1:length(ν)
            name_ν  = @sprintf("%.3fMHz",ν[β]/1e6)
            if !(name_ν in names(file))
                group_ν = g_create(file,name_ν)
            else
                group_ν = file[name_ν]
            end

            group_m = g_create(group_ν,string(m))
            group_m["block"] = v[β].block
        end
    end
end

function load_mmodes(filename, ν)
    blocks = MModeBlock[]
    jldopen(filename,"r") do file
        mmax = file["mmax"] |> read

        name_ν  = @sprintf("%.3fMHz",ν/1e6)
        group_ν = file[name_ν]

        for m = 0:mmax
            group_m = group_ν[string(m)]
            block = group_m["block"] |> read
            push!(blocks,MModeBlock(block,m,ν))
        end
    end
    MModes(blocks)
end

function load_mmodes(filename, ν, m::Int)
    blocks = MModeBlock[]
    jldopen(filename,"r") do file
        mmax = file["mmax"] |> read

        for β = 1:length(ν)
            name_ν  = @sprintf("%.3fMHz",ν[β]/1e6)
            group_ν = file[name_ν]

            group_m = group_ν[string(m)]
            block = group_m["block"] |> read
            push!(blocks,MModeBlock(block,m,ν[β]))
        end
    end
    MModes(blocks)
end

