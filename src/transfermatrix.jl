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

immutable TransferMatrix{Nbase,lmax,mmax}
    blocks::Vector{Matrix{Complex128}}
end

function TransferMatrix(Nbase,lmax,mmax)
    blocks = [zeros(Complex128,two(m)*Nbase,lmax-m+1) for m = 0:mmax]
    TransferMatrix{Nbase,lmax,mmax}(blocks)
end

blocks(B::TransferMatrix) = B.blocks
Nbase{Nbase,lmax,mmax}(B::TransferMatrix{Nbase,lmax,mmax}) = Nbase
lmax{Nbase,lmax,mmax}(B::TransferMatrix{Nbase,lmax,mmax}) = lmax
mmax{Nbase,lmax,mmax}(B::TransferMatrix{Nbase,lmax,mmax}) = mmax

# α labels the baseline
# l and m label the spherical harmonic

block(B::TransferMatrix,m) = blocks(B)[abs(m)+1]
getindex(B::TransferMatrix,α,l,m) = block(B,m)[α+one(m)*Nbase(B),l-abs(m)+1]
setindex!(B::TransferMatrix,x,α,l,m) = block(B,m)[α+one(m)*Nbase(B),l-abs(m)+1] = x

function TransferMatrix(beam::HEALPixMap,
                        positions::Matrix{Float64},
                        frequency::Float64,
                        phasecenter::NTuple{3,Float64};
                        lmax::Int = 100,
                        mmax::Int = 100)
    Nant  = size(positions,2)
    Nbase = div(Nant*(Nant-1),2)
    B = TransferMatrix(Nbase,lmax,mmax)
    populate!(B,beam,positions,frequency,phasecenter)
    B
end

function populate!(B::TransferMatrix,
                   beam, positions, frequency, phasecenter)
    Nant = round(Int,(1+sqrt(1+8Nbase(B)))/2)
    λ = c / frequency
    α = 1
    for ant1 = 1:Nant, ant2 = ant1+1:Nant
        @show ant1,ant2
        u = (positions[1,ant1] - positions[1,ant2])/λ
        v = (positions[2,ant1] - positions[2,ant2])/λ
        w = (positions[3,ant1] - positions[3,ant2])/λ
        # Account for the extra phase due to the fact that the phase
        # center is likely not perfectly orthogonal to the baseline.
        extraphase = -2π*(u*phasecenter[1]+v*phasecenter[2]+w*phasecenter[3])
        # Use spherical harmonic transforms to incorporate the effects
        # of the beam.
        realfringe,imagfringe = planewave(u,v,w,extraphase,lmax=lmax(B),mmax=mmax(B))
        realbeamfringe = map2alm(beam.*alm2map(realfringe,nside=512),lmax=lmax(B),mmax=mmax(B))
        imagbeamfringe = map2alm(beam.*alm2map(imagfringe,nside=512),lmax=lmax(B),mmax=mmax(B))
        # Pack the transfer matrix
        # (the conjugations come about because Shaw et al. 2014, 2015
        # actually expand the baseline pattern in terms of the
        # spherical harmonic conjugates)
        for l = 0:lmax(B)
            B[α,l,0] = conj(realbeamfringe[l,0]) + 1im*conj(imagbeamfringe[l,0])
        end
        for m = 1:mmax(B), l = m:lmax(B)
            B[α,l,+m] = conj(realbeamfringe[l,m]) + 1im*conj(imagbeamfringe[l,m])
            B[α,l,-m] = conj(realbeamfringe[l,m]) - 1im*conj(imagbeamfringe[l,m])
        end
        α += 1
    end
    nothing
end

