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

one(m) = ifelse(m  > 0, 1, 0)
two(m) = ifelse(m != 0, 2, 1)

immutable TransferMatrix{Nbase,Nfreq,lmax,mmax}
    blocks::Matrix{Matrix{Complex128}}
end

function TransferMatrix(Nbase,Nfreq,lmax,mmax)
    blocks = [zeros(Complex128,two(m)*Nbase,lmax-m+1) for β = 1:Nfreq, m = 0:mmax]
    TransferMatrix{Nbase,Nfreq,lmax,mmax}(blocks)
end

blocks(B::TransferMatrix) = B.blocks
Nbase{Nbase,Nfreq,lmax,mmax}(B::TransferMatrix{Nbase,Nfreq,lmax,mmax}) = Nbase
Nfreq{Nbase,Nfreq,lmax,mmax}(B::TransferMatrix{Nbase,Nfreq,lmax,mmax}) = Nfreq
lmax{Nbase,Nfreq,lmax,mmax}(B::TransferMatrix{Nbase,Nfreq,lmax,mmax}) = lmax
mmax{Nbase,Nfreq,lmax,mmax}(B::TransferMatrix{Nbase,Nfreq,lmax,mmax}) = mmax

# α labels the baseline
# β labels the frequency
# l and m label the spherical harmonic

function getindex(B::TransferMatrix,α,β,l,m)
    blocks(B)[β,abs(m)+1][α+one(m)*Nbase(B),l-abs(m)+1]
end

function setindex!(B::TransferMatrix,x,α,β,l,m)
    blocks(B)[β,abs(m)+1][α+one(m)*Nbase(B),l-abs(m)+1] = x
end

function TransferMatrix(beam::Alm,
                        positions::Matrix{Float64},
                        frequencies::Vector{Float64},
                        cg::CGTable;
                        lmax::Int = 100,
                        mmax::Int = 100)
    Nant  = size(positions,2)
    Nbase = div(Nant*(Nant-1),2)
    Nfreq = length(frequencies)
    B = TransferMatrix(Nbase,Nfreq,lmax,mmax)
    populate!(B,beam,positions,frequencies,cg)
    B
end

function populate!(B::TransferMatrix,
                   positions, frequencies, cg)
    for β = 1:Nfreq
        ν = frequencies[β]
        λ = TTCal.c / ν
        α = 1
        for ant1 = 1:Nant, ant2 = ant1+1:Nant
            u = (antenna_positions[1,ant2] - antenna_positions[1,ant1])/λ
            v = (antenna_positions[2,ant2] - antenna_positions[2,ant1])/λ
            w = (antenna_positions[3,ant2] - antenna_positions[3,ant1])/λ
            realfringe,imagfringe = planewave(u,v,w,lmax=lmax,mmax=mmax)
            realbeamfringe = multiply(beam,realfringe)
            imagbeamfringe = multiply(beam,imagfringe)
            # Pack the transfer matrix
            # (the conjugations come about because Shaw et al. 2014, 2015
            # actually expands the baseline pattern in terms of the
            # spherical harmonic conjugates)
            for l = 0:lmax
                B[α,β,l,0] = conj(realbeamfringe[l,0]) + 1im*conj(imagebeamfringe[l,0])
            end
            for m = 1:mmax, l = m:lmax
                B[α,β,l,+m] = conj(realbeamfringe[l,m]) + 1im*conj(imagebeamfringe[l,m])
                B[α,β,l,-m] = conj(realbeamfringe[l,m]) - 1im*conj(imagebeamfringe[l,m])
            end
            α += 1
        end
    end
    nothing
end

