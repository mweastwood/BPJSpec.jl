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

Nfreq(C::CovarianceMatrix) = size(C.blocks,1)
blocksize(C::CovarianceMatrix) = size(C[1,1].block,1)
getindex(C::CovarianceMatrix,β1,β2) = C.blocks[β1,β2]

"""
    congruence(B::SpectralTransferMatrix,C::CovarianceMatrix) -> B*C*B'

Compute the congruence transform of the covariance matrix with respect to
the transfer matrix. That is, change the basis of the covariance matrix
from spherical harmonic coefficients to m-modes.
"""
function congruence(B::SpectralTransferMatrix,C::CovarianceMatrix)
    # TODO: write this
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
    (A * ((l+1)/100)^(-α)
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

