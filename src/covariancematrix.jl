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
    CovarianceMatrixBlock

A single block of a `CovarianceMatrix` corresponding to a pair
of frequency channels.
"""
immutable CovarianceMatrixBlock
    block::Matrix{Complex128}
end

"""
    CovarianceMatrix

A covariance matrix where frequency channels are arranged in arranged
in blocks.
"""
immutable CovarianceMatrix
    blocks::Matrix{CovarianceMatrixBlock}
end

function CovarianceMatrix(Nfreq)
    CovarianceMatrix(Array{CovarianceMatrixBlock}(Nfreq,Nfreq))
end

Nfreq(C::CovarianceMatrix) = size(C.blocks,1)
blocksize(C::CovarianceMatrix) = size(C[1,1].block,1)
getindex(C::CovarianceMatrix,β1,β2) = C.blocks[β1,β2]
setindex!(C::CovarianceMatrix,x,β1,β2) = C.blocks[β1,β2] = x

"""
    congruence(B::SpectralTransferMatrix,C::CovarianceMatrix) -> B*C*B'

Compute the congruence transform of the covariance matrix with respect to
the transfer matrix. That is, change the basis of the covariance matrix
from spherical harmonic coefficients to m-modes.
"""
function congruence(B::SpectralTransferMatrix,C::CovarianceMatrix)
    Nfreq(B) == Nfreq(C) || error("The values of Nfreq must be the same.")
    out = CovarianceMatrix(Nfreq(C))
    for β2 = 1:Nfreq(C), β1 = 1:Nfreq(C)
        out[β1,β2] = CovarianceMatrixBlock(B[β1].block*C[β1,β2].block*B[β2].block')
    end
    out
end

abstract AbstractComponent

"""
    CovarianceMatrix(component::AbstractComponent,ν,lmax,m)

Construct a covariance matrix for the given component of the sky.
"""
function CovarianceMatrix(component::AbstractComponent,ν,lmax,m)
    Nfreq = length(ν)
    blocksize = lmax-m+1
    blocks = Array{CovarianceMatrixBlock}(Nfreq,Nfreq)
    for β2 = 1:Nfreq, β1 = 1:Nfreq
        block = zeros(Complex128,blocksize,blocksize)
        for l = m:lmax
            block[l-m+1,l-m+1] = component(l,ν[β1],ν[β2])
        end
        blocks[β1,β2] = CovarianceMatrixBlock(block)
    end
    CovarianceMatrix(blocks)
end

function Base.full(C::CovarianceMatrix)
    N = blocksize(C)*Nfreq(C)
    fullC = zeros(Complex128,N,N)
    for β2 = 1:Nfreq(C), β1 = 1:Nfreq(C)
        idx1 = (β1-1)*blocksize(C)+1:β1*blocksize(C)
        idx2 = (β2-1)*blocksize(C)+1:β2*blocksize(C)
        subC = sub(fullC,idx1,idx2)
        block = C[β1,β2].block
        subC[:] = block
    end
    fullC
end

"""
    Cforeground(l,ν1,ν2,ν0,A,α,β,ζ) -> Cl::Float64

Evaluate a model for the multifrequency angular power spectrum of
a foreground component. Note that `ν0`, `A`, `α`, `β`, and `ζ` are
all parameters of the model.
"""
function Cforeground(l,ν1,ν2,ν0,A,α,β,ζ)
    (A * (l+1)^(-α)
        * (ν1*ν2/ν0^2)^(-β)
        * exp(-log(ν1/ν2)^2/(2*ζ^2)))
end

immutable ForegroundModel <: AbstractComponent
    ν0::Float64
    A::Float64
    α::Float64
    β::Float64
    ζ::Float64
end

function call(model::ForegroundModel,l,ν1,ν2)
    Cforeground(l,ν1,ν2,model.ν0,model.A,model.α,model.β,model.ζ)
end

"""
    Csignal_spherical(l,ν1,ν2,kedges,P) -> Cl::Float64

Project a model of a spherically averaged spatial power spectrum
to calculate the corresponding multifrequency angular power spectrum.
There are two branches because a limit must be taken as `ν1-ν2` tends
to zero.

This projection assumes that the power spectrum is constant within
the bins defined by `kedges`. The integral is then analytically
evaluated within each bin.

    1/(π*r1*r2) * ∫dk*P(k)*cos(sqrt(k^2-k⟂^2)*Δr)*k/sqrt(k^2-k⟂^2)

Note that because `kedges` defines the edges of each of the bins,
the length of `kedges` should be one greater than the length of `P`.
"""
function Csignal_spherical(l,ν1,ν2,kedges,P)
    if ν1 == ν2
        z = redshift(ν1)
        r = comoving_distance(z)
        kperp = l/r
        out = 0.0
        for i = 1:length(P)
            kperp ≥ kedges[i+1] && continue
            kmax = kedges[i+1]
            kmin = max(kedges[i],kperp)
            out += (sqrt(kmax^2-kperp^2)
                   -sqrt(kmin^2-kperp^2))*P[i]/(π*r^2)
        end
        return out
    else
        z1 = redshift(ν1)
        z2 = redshift(ν2)
        r1 = comoving_distance(z1)
        r2 = comoving_distance(z2)
        Δr = r2-r1
        rmean = (r1+r2)/2
        kperp = l/rmean
        out = 0.0
        for i = 1:length(P)
            kperp ≥ kedges[i+1] && continue
            kmax = kedges[i+1]
            kmin = max(kedges[i],kperp)
            out += (sin(sqrt(kmax^2-kperp^2)*Δr)
                   -sin(sqrt(kmin^2-kperp^2)*Δr))*P[i]/(π*Δr*r1*r2)
        end
        return out
    end
end

immutable SphericalSignalModel <: AbstractComponent
    kedges::Vector{Float64} # These define the edges of power spectrum bins
    power::Vector{Float64}
end

function call(model::SphericalSignalModel,l,ν1,ν2)
    Csignal_spherical(l,ν1,ν2,model.kedges,model.power)
end

# TODO: implement CylindricalSignalModel
# This will have cylindrical binning of the power spectrum.
# That is, allow the power spectrum to be specified as a
# function of kperp and kpara.

