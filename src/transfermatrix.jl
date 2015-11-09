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
    healpix(frame::ReferenceFrame, beam::BeamModel, frequency)

Create a Healpix map of a TTCal beam model.

At the moment only the $I \rightarrow I$ element of the Mueller
matrix is used.
"""
function healpix(frame::ReferenceFrame, beam::BeamModel, frequency)
    zenith = measure(frame,Direction(dir"AZEL",Quantity(0.0,"deg"),Quantity(90.0,"deg")),dir"APP")
    zenith_dec = latitude(zenith,"rad")
    zenith_vec = LibHealpix.ang2vec(π/2-zenith_dec,0.0)
    global_north = [0.0,0.0,1.0]
    local_north  = gramschmidt(global_north,zenith_vec)
    map = HealpixMap(zeros(nside2npix(512)))
    for i = 1:length(map)
        vec = LibHealpix.pix2vec_ring(512,i)
        el  = π/2 - angle_between(vec,zenith_vec)
        vec = gramschmidt(vec,zenith_vec)
        az  = angle_between(vec,local_north)
        if el < 0
            map[i] = 0
        else
            J = beam(frequency,az,el)
            M = mueller(J)
            map[i] = M[1,1]
        end
    end
    map
end

doc"""
    immutable TransferMatrixBlock <: MatrixBlock

This type stores a single block of the transfer matrix.

**Fields:**

* `block` stores the actual matrix block
* `lmax` stores the maximum value of $l$
* `m` stores the value of $m$ corresponding to this block of the matrix
* `ν` stores the frequency (in Hz)
"""
immutable TransferMatrixBlock <: MatrixBlock
    block::Matrix{Complex128}
    lmax::Int
    m::Int
    ν::Float64
end

default_size(::Type{TransferMatrixBlock},Nbase,lmax,m) = (two(m)*Nbase, lmax-m+1)

"""
    TransferMatrixBlock(Nbase,lmax,m,ν)

Create an empty transfer matrix block of the default size.
"""
function TransferMatrixBlock(Nbase,lmax,m,ν)
    sz = default_size(TransferMatrixBlock,Nbase,lmax,m)
    block = zeros(Complex128,sz)
    TransferMatrixBlock(block,lmax,m,ν)
end

@enum Representation one_ν one_m

doc"""
    immutable TransferMatrix{sym}

The transfer matrix represents the instrumental response of the
interferometer to the spherical harmonic coefficients of the sky.

The transfer matrix is usually very large, but we can use its block
diagonal structure to work with parts of the matrix separately.
We can choose between:

1. Working with one frequency channel at a time, but all values of $m$,
   where $m$ is the azimuthal quantum number of the spherical harmonics, or
2. Working with several frequency channels, but only one value of $m$.

The former representation is more useful while calculating elements of
the transfer matrix, but the latter representation is more useful for
foreground filtering and cosmological power spectrum estimation (where
frequency structure is important). The type parameter `sym` is therefore
the switch that indicates whether we are working with the first
representation or the second.

    TransferMatrix{one_ν} # one frequency channel, many m
    TransferMatrix{one_m} # many frequency channels, one m

**Fields:**

* `blocks` is a list of `TransferMatrixBlock`s
* `lmax` and `mmax` stores the maximum values of $l$ and $m$ (the spherical
   harmonic quantum numbers) respectively.
"""
immutable TransferMatrix{rep}
    blocks::Vector{TransferMatrixBlock}
end

"""
    TransferMatrix(blocks)

Construct a transfer matrix from the list of blocks.

The representation (`one_ν` or `one_m`) is inferred from the provided
list of blocks.
"""
function TransferMatrix(blocks)
    lmax = [block.lmax for block in blocks]
    length(unique(lmax)) == 1 || error("All blocks in a TransferMatrix must have a consistent lmax.")

    m = [block.m for block in blocks]
    ν = [block.ν for block in blocks]
    mmax = length(m) - 1

    # if all the blocks are the same frequency and are sorted in order
    # of increasing m, we have the one_ν representation of the transfer
    # matrix
    if length(unique(ν)) == 1 && m == collect(0:mmax)
        return TransferMatrix{one_ν}(blocks)
    end

    # if all the blocks have different frequencies (not necessarily sorted)
    # but only have one value of m, we have the one_m representation of
    # the transfer matrix
    if length(unique(ν)) == length(ν) && length(unique(m)) == 1
        return TransferMatrix{one_m}(blocks)
    end

    error("""
    The list of blocks must either:

    * have a consistent frequency channel representing all values of m from 0 to some mmax, or
    * have different frequency channels, but a consistent value of m
    """)
end

# α labels the baseline
# β labels the frequency channel
# l and m label the spherical harmonic

getindex(B::TransferMatrixBlock,α,l) = B.block[α,l-B.m+1]
setindex!(B::TransferMatrixBlock,x,α,l) = B.block[α,l-B.m+1] = x

getindex(B::TransferMatrix{one_ν},m) = B.blocks[m+1]
getindex(B::TransferMatrix{one_ν},α,l,m) = B[m][α,l]
setindex!(B::TransferMatrix{one_ν},x,α,l,m) = B[m][α,l] = x

getindex(B::TransferMatrix{one_m},β) = B.blocks[β]
getindex(B::TransferMatrix{one_m},α,β,l) = B[β][α,l]
setindex!(B::TransferMatrix{one_m},x,α,β,l) = B[β][α,l] = x

#=
function TransferMatrix(Nbase,lmax,mmax)
    blocks = [TransferMatrixBlock(Nbase,lmax,mmax,m) for m = 0:mmax]
    TransferMatrix{lmax,mmax}(blocks)
end


function SpectralTransferMatrix(Nbase,Nfreq,lmax,mmax,m)
    blocks = [TransferMatrixBlock(Nbase,lmax,mmax,m) for β = 1:Nfreq]
    SpectralTransferMatrix{lmax,mmax}(m,blocks)
end

for T in (:TransferMatrixBlock,:TransferMatrix,:SpectralTransferMatrix)
    @eval lmax{l,m}(B::$T{l,m}) = l
    @eval mmax{l,m}(B::$T{l,m}) = m
end

Nfreq(B::SpectralTransferMatrix) = length(B.blocks)

function ==(lhs::TransferMatrixBlock,rhs::TransferMatrixBlock)
    lhs.block == rhs.block && lhs.m == rhs.m &&
        lmax(lhs) == lmax(rhs) && mmax(lhs) == mmax(rhs)
end

function ==(lhs::TransferMatrix,rhs::TransferMatrix)
    lhs.blocks == rhs.blocks && lmax(lhs) == lmax(rhs) && mmax(lhs) == mmax(rhs)
end

Base.svd(B::TransferMatrixBlock) = svd(B.block)

Base.size(B::TransferMatrixBlock) = size(B.block)
function Base.size(B::TransferMatrix)
    x = 0; y = 0
    for m = 0:mmax(B)
        sz = size(B[m])
        x += sz[1]
        y += sz[2]
    end
    x,y
end
function Base.size(B::SpectralTransferMatrix)
    x = 0; y = 0
    for β = 1:Nfreq(B)
        sz = size(B[β])
        x += sz[1]
        y += sz[2]
    end
    x,y
end
=#

"""
    gentransfer(beam::HealpixMap, u, v, w, ν, phasecenter; lmax = 100, mmax = 100)

Construct a transfer matrix for the given beam model and observational parameters.
"""
function gentransfer(beam::HealpixMap, u, v, w, ν, phasecenter;
                     lmax::Int = 100, mmax::Int = 100)
    Nbase = length(u)
    B = [TransferMatrixBlock(Nbase,lmax,m,ν) for m = 0:mmax] |> TransferMatrix
    gentransfer!(B,beam,u,v,w,ν,phasecenter,lmax=lmax,mmax=mmax)
    B
end

function gentransfer!(B::TransferMatrix{one_ν},
                      beam::HealpixMap, u, v, w, ν, phasecenter;
                      lmax::Int = 100, mmax::Int = 100)
    Nbase = length(u)
    λ = c / ν
    @showprogress 1 "Computing transfer matrix..." for α = 1:Nbase
        uλ = u[α]/λ
        vλ = v[α]/λ
        wλ = w[α]/λ
        # Account for the extra phase due to the fact that the phase
        # center is likely not perfectly orthogonal to the baseline.
        extraphase = -2π*(uλ*phasecenter[1]+vλ*phasecenter[2]+wλ*phasecenter[3])
        realfringe,imagfringe = planewave(uλ,vλ,wλ,extraphase,lmax=lmax,mmax=mmax)
        # Use spherical harmonic transforms to incorporate the effects of the beam.
        realbeamfringe = map2alm(beam.*alm2map(realfringe,nside=512),lmax=lmax,mmax=mmax)
        imagbeamfringe = map2alm(beam.*alm2map(imagfringe,nside=512),lmax=lmax,mmax=mmax)
        # Pack the transfer matrix
        # (the conjugations come about because Shaw et al. 2014, 2015 actually expand the
        # baseline pattern in terms of the spherical harmonic conjugates)
        for l = 0:lmax
            B[α,l,0] = conj(realbeamfringe[l,0]) + 1im*conj(imagbeamfringe[l,0])
        end
        for m = 1:mmax, l = m:lmax
            α1 = α         # positive m
            α2 = α + Nbase # negative m
            B[α1,l,m] = conj(realbeamfringe[l,m]) + 1im*conj(imagbeamfringe[l,m])
            B[α2,l,m] = conj(realbeamfringe[l,m]) - 1im*conj(imagbeamfringe[l,m])
        end
        α += 1
    end
    B
end

#=
"""
    SpectralTransferMatrix(m,matrices::Vector{TransferMatrix})

Construct a SpectralTransferMatrix from the given list of transfer matrices.
Each transfer matrix should correspond to a different frequency channel.
"""
function SpectralTransferMatrix(m,matrices::Vector{TransferMatrix})
    lmax′ = lmax(matrices[1])
    mmax′ = mmax(matrices[1])
    blocks = TransferMatrixBlock[]
    for matrix in matrices
        if lmax(matrix) != lmax′ || mmax(matrix) != mmax′
            error("The transfer matrices must all have the same lmax and mmax.")
        end
        push!(blocks,matrix[m])
    end
    SpectralTransferMatrix{lmax′,mmax′}(m,blocks)
end

function Base.full(B::SpectralTransferMatrix)
    out = zeros(Complex128,size(B))
    idx1 = 1; idx2 = 1
    for β = 1:Nfreq(B)
        block = B[β].block
        out[idx1:idx1+size(block,1)-1,
            idx2:idx2+size(block,2)-1] = block
        idx1 += size(block,1)
        idx2 += size(block,2)
    end
    out
end
=#

