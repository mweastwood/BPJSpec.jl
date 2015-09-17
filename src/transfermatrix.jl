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

"""
    TransferMatrixBlock{lmax,mmax}

This type stores a single block of the transfer matrix. These blocks
correspond to a single frequency channel and a single value of `m`.
"""
immutable TransferMatrixBlock{lmax,mmax}
    m::Int
    block::Matrix{Complex128}
end

function TransferMatrixBlock(Nbase,lmax,mmax,m)
    TransferMatrixBlock{lmax,mmax}(m,zeros(Complex128,two(m)*Nbase,lmax-m+1))
end

"""
    TransferMatrix{lmax,max}

This type stores all blocks of the transfer matrix corresponding to
a single frequency channel, but all values of `m`.
"""
immutable TransferMatrix{lmax,mmax}
    blocks::Vector{TransferMatrixBlock{lmax,mmax}}
end

function TransferMatrix(Nbase,lmax,mmax)
    blocks = [TransferMatrixBlock(Nbase,lmax,mmax,m) for m = 0:mmax]
    TransferMatrix{lmax,mmax}(blocks)
end

"""
    SpectralTransferMatrix{lmax,max}

This type stores all blocks of the transfer matrix corresponding to
a single value of `m`, but all frequency channels.
"""
immutable SpectralTransferMatrix{lmax,mmax}
    m::Int
    blocks::Vector{TransferMatrixBlock{lmax,mmax}}
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

# α labels the baseline
# β labels the frequency channel
# l and m label the spherical harmonic

getindex(B::TransferMatrixBlock,α,l) = B.block[α,l-B.m+1]
setindex!(B::TransferMatrixBlock,x,α,l) = B.block[α,l-B.m+1] = x

getindex(B::TransferMatrix,m) = B.blocks[m+1]
getindex(B::TransferMatrix,α,l,m) = B[m][α,l]
setindex!(B::TransferMatrix,x,α,l,m) = B[m][α,l] = x

getindex(B::SpectralTransferMatrix,β) = B.blocks[β]
getindex(B::SpectralTransferMatrix,α,β,l) = B[β][α,l]
setindex!(B::SpectralTransferMatrix,x,α,β,l) = B[β][α,l] = x

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

"""
    TransferMatrix(beam::HealpixMap,obs::ObsParam,channel;lmax=100,mmax=100)

Construct a transfer matrix for the given beam model and observational parameters.
"""
function TransferMatrix(beam::HealpixMap,
                        obs::ObsParam,
                        channel;
                        lmax::Int = 100,
                        mmax::Int = 100)
    positions = antpos(obs)
    frequency = freq(obs)[channel]
    phase     = phasecenter(obs)
    B = TransferMatrix(Nbase(obs),lmax,mmax)
    populate!(B,beam,obs,channel)
    B
end

function populate!(B::TransferMatrix,
                   beam::HealpixMap,
                   obs::ObsParam,
                   channel)
    positions = antpos(obs)
    λ = c / freq(obs)[channel]
    phase = phasecenter(obs)
    α = 1
    for ant1 = 1:Nant(obs), ant2 = ant1+1:Nant(obs)
        @show ant1,ant2
        u = (positions[1,ant1] - positions[1,ant2])/λ
        v = (positions[2,ant1] - positions[2,ant2])/λ
        w = (positions[3,ant1] - positions[3,ant2])/λ
        # Account for the extra phase due to the fact that the phase
        # center is likely not perfectly orthogonal to the baseline.
        extraphase = -2π*(u*phase[1]+v*phase[2]+w*phase[3])
        # Use spherical harmonic transforms to incorporate the effects
        # of the beam.
        realfringe,imagfringe = planewave(u,v,w,extraphase,lmax=lmax(B),mmax=mmax(B))
        realbeamfringe = map2alm(beam.*alm2map(realfringe,nside=512),
                                 lmax=lmax(B),mmax=mmax(B))
        imagbeamfringe = map2alm(beam.*alm2map(imagfringe,nside=512),
                                 lmax=lmax(B),mmax=mmax(B))
        # Pack the transfer matrix
        # (the conjugations come about because Shaw et al. 2014, 2015
        # actually expand the baseline pattern in terms of the
        # spherical harmonic conjugates)
        for l = 0:lmax(B)
            B[α,l,0] = conj(realbeamfringe[l,0]) + 1im*conj(imagbeamfringe[l,0])
        end
        for m = 1:mmax(B), l = m:lmax(B)
            α1 = α              # positive m
            α2 = α + Nbase(obs) # negative m
            B[α1,l,m] = conj(realbeamfringe[l,m]) + 1im*conj(imagbeamfringe[l,m])
            B[α2,l,m] = conj(realbeamfringe[l,m]) - 1im*conj(imagbeamfringe[l,m])
        end
        α += 1
    end
    nothing
end

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

