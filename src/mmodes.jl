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
    MModeBlock{mmax}

This type stores the m-modes corresponding to a single block of
the transfer matrix.
"""
immutable MModeBlock{mmax}
    m::Int
    block::Vector{Complex128}
end

function MModeBlock(Nbase,mmax,m)
    MModeBlock{mmax}(m,zeros(Complex128,two(m)*Nbase))
end

"""
    MModes{mmax}

This type stores the m-modes corresponding to a single frequency
channel, but all values of `m`.
"""
immutable MModes{mmax}
    blocks::Vector{MModeBlock{mmax}}
end

function MModes(Nbase,mmax)
    blocks = [MModeBlock(Nbase,mmax,m) for m = 0:mmax]
    MModes{mmax}(blocks)
end

"""
    SpectralMModes{mmax}

This type stores the m-modes corresponding to a single value of
`m`, but all frequency channels.
"""
immutable SpectralMModes{mmax}
    blocks::Vector{MModeBlock{mmax}}
end

function SpectralMModes(Nbase,Nfreq,mmax,m)
    blocks = [MModeBlock(Nbase,mmax,m) for β = 1:Nfreq]
    MModes{mmax}(blocks)
end

for T in (:MModeBlock,:MModes,:SpectralMModes)
    @eval mmax{m}(v::$T{m}) = m
end

function ==(lhs::MModeBlock,rhs::MModeBlock)
    lhs.block == rhs.block && lhs.m == rhs.m && mmax(lhs) == mmax(rhs)
end

function ==(lhs::MModes,rhs::MModes)
    lhs.blocks == rhs.blocks && mmax(lhs) == mmax(rhs)
end

# α labels the baseline
# β labels the frequency channel
# m labels the mode

getindex(v::MModeBlock,α) = v.block[α]
setindex!(v::MModeBlock,x,α) = v.block[α] = x

getindex(v::MModes,m) = v.blocks[m+1]
getindex(v::MModes,α,m) = v[m][α]
setindex!(v::MModes,x,α,m) = v[m][α] = x

getindex(v::SpectralMModes,β) = v.blocks[β]
getindex(v::SpectralMModes,α,β) = v.blocks[β][α]
setindex!(v::SpectralMModes,x,α,β) = v.blocks[β][α] = x

function *(B::TransferMatrix,alm::Alm)
    lmax(B) == lmax(alm) || error("The values of lmax must be the same.")
    mmax(B) == mmax(alm) || error("The values of mmax must be the same.")
    blocks = [MModeBlock{mmax(B)}(m,B[m].block*block(alm,m)) for m = 0:mmax(B)]
    MModes{mmax(B)}(blocks)
end

"""
    MModes(visibilities;mmax=100)

Calculate the m-modes from the given visibilities. The visibilities
should be provided as a matrix where the first dimension indicates
the baseline and the second dimension indicates the sidereal time.
That is, the FFT will be performed over the second dimension. The
visibilities should span a full sidereal day.
"""
function MModes{T<:Complex}(visibilities::Matrix{T};
                            mmax::Int = 100)
    Nbase,Ntime = size(visibilities)
    M = fft(visibilities,2)/Ntime
    v = MModes(Nbase,mmax)
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
    visibilities(v::MModes)

Calculate the visibilities from the given m-modes. The visibilities
will be returned as a two dimensional array where the first dimension
indicates the baseline and the second dimension indicates the sidereal
time, which evenly spans a full sidereal day.
"""
function visibilities(v::MModes)
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

