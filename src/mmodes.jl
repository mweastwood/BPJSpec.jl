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

immutable MModes{Nbase,mmax}
    blocks::Vector{Vector{Complex128}}
end

function MModes(Nbase,mmax)
    blocks = [zeros(Complex128,two(m)*Nbase) for m = 0:mmax]
    MModes{Nbase,mmax}(blocks)
end

blocks(v::MModes) = v.blocks
Nbase{Nbase,mmax}(v::MModes{Nbase,mmax}) = Nbase
mmax{Nbase,mmax}(v::MModes{Nbase,mmax}) = mmax

function ==(lhs::MModes,rhs::MModes)
    blocks(lhs) == blocks(rhs) &&
        Nbase(lhs) == Nbase(rhs) &&
        mmax(lhs) == mmax(rhs)
end

# α labels the baseline
# m labels the mode

block(v::MModes,m) = blocks(v)[abs(m)+1]
getindex(v::MModes,α,m) = block(v,m)[α+one(m)*Nbase(v)]
setindex!(v::MModes,x,α,m) = block(v,m)[α+one(m)*Nbase(v)] = x

function MModes(transfermatrix::TransferMatrix,alm::Alm)
    mmodes = MModes(Nbase(transfermatrix),
                    mmax(transfermatrix))
    for m = 0:mmax(transfermatrix)
        B = block(transfermatrix,m)
        a = block(alm,m)
        v = B*a
        mmodes.blocks[m+1] = v
    end
    mmodes
end

*(transfermatrix::TransferMatrix,alm::Alm) = MModes(transfermatrix,alm)

function MModes{T<:Complex}(visibilities::Matrix{T};mmax::Int = 100)
    Nbase,Ntime = size(visibilities)
    M = fft(visibilities,2)/Ntime

    mmodes = MModes(Nbase,mmax)
    for α = 1:Nbase
        mmodes[α,0] = M[α,1]
    end
    for m = 1:mmax, α = 1:Nbase
        mmodes[α,+m] =      M[α,m+1]
        mmodes[α,-m] = conj(M[α,Ntime+1-m])
    end
    mmodes
end

function visibilities(mmodes::MModes)
    M = zeros(Complex128,Nbase(mmodes),2mmax(mmodes)+1)
    for α = 1:Nbase(mmodes)
        M[α,1] = mmodes[α,0]
    end
    for m = 1:mmax(mmodes), α = 1:Nbase(mmodes)
        M[α,m+1]               =      mmodes[α,+m]
        M[α,2mmax(mmodes)+2-m] = conj(mmodes[α,-m])
    end
    ifft(M,2)*size(M,2)
end

