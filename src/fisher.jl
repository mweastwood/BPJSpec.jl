# Copyright (c) 2015-2017 Michael Eastwood
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

struct FisherMatrix
    matrix :: BlockDiagonalMatrix
end

# The disk I/O is killing me. It's so slow, but the angular covariance matrices are small enough to
# fit in memory. So we'll use this type to store the pre-loaded blocks.
struct BasisCovarianceMatrices
    lmax     :: Int
    matrices :: Vector{AngularCovarianceMatrix}
    blocks   :: Matrix{Matrix{Float64}} # pre-loaded blocks from all the matrices
end

function preload(matrices::Vector{AngularCovarianceMatrix})
    N = length(matrices)
    lmax = matrices[1].lmax
    blocks = Array{Matrix{Float64}}(lmax+1, N)
    for idx = 1:N, l = 0:lmax
        blocks[l+1, idx] = matrices[idx][l]
    end
    BasisCovarianceMatrices(lmax, matrices, blocks)
end

function FisherMatrix(transfermatrix, covariancematrices)
    #@time basis = preload(covariancematrices)
    #save("/dev/shm/mweastwood/cache.jld2", "basis", basis)
    @time basis = load("/dev/shm/mweastwood/cache.jld2", "basis")
    compute_fisher_matrix(transfermatrix, basis, 100)
end

function compute_fisher_matrix(transfermatrix, basis, m)
    N  = length(basis.matrices)
    F  = zeros(N, N) # the Fisher matrix
    @time BB = reorder(transfermatrix, m)
    BBCa = [similar(BB[idx, jdx]) for jdx = 1:size(BB, 2), idx = 1:size(BB, 1)]
    BBCb = [similar(BB[idx, jdx]) for jdx = 1:size(BB, 2), idx = 1:size(BB, 1)]

    prg = Progress((N*(N+1))÷2)
    for a = 1:N
        Ca = view(basis.blocks, m+1:basis.lmax+1, a)
        multiply!(BBCa, BB, Ca)
        F[a, a] = real(tr(BBCa, BBCa))
        next!(prg)

        #for b = a+1:N
        #    Cb = view(basis.blocks, m+1:basis.lmax+1, b)
        #    multiply!(BBCb, BB, Cb)
        #    F[a, b] = real(tr(BBCa, BBCb))
        #    F[b, a] = F[a, b]
        #    next!(prg)
        #end
    end
    F
end

function reorder(transfermatrix::SpectralBlockDiagonalMatrix, m)
    # The angular covariance matrix is block diagonal in l. The transfer matrix is not, but we can
    # still probably save some computer time by splitting the transfer matrix into its corresponding
    # blocks. This will hopefully reduce the amount of time we spend operating on the structural
    # zeroes of the covariance matrix.
    lmax  = transfermatrix.mmax
    Nfreq = length(transfermatrix.frequencies)

    blocks = [transfermatrix[m, ν] for ν in transfermatrix.frequencies]
    blocks = [block'*block for block in blocks]

    N = lmax - m + 1
    output = Array{Matrix{Complex128}}(N, N)
    for l1 = m:lmax, l2 = m:lmax
        output_block = zeros(Complex128, Nfreq, Nfreq)
        for β = 1:Nfreq
            block = blocks[β]
            output_block[β, β] = block[l1-m+1, l2-m+1]
        end
        output[l1-m+1, l2-m+1] = output_block
    end
    output
end

function reorder(covariancematrix::AngularCovarianceMatrix, m)
    lmax  = covariancematrix.lmax
    Nfreq = length(covariancematrix.frequencies)
    [covariancematrix[l] for l = m:lmax]
end

function multiply!(output, BB, C)
    N = length(C)
    for jdx = 1:N, idx = 1:N
        A_mul_B!(output[idx, jdx], BB[idx, jdx], C[idx])
    end
    output
end

"""
    tr(A, B) -> tr(AB)

The trace of the product of two matrices, computed without first computing the product itself.
"""
function tr(A, B)
    trace = zero(eltype(A))
    @inbounds for j = 1:size(A, 1), i = 1:size(A, 2)
        trace += A[i, j] * B[j, i]
    end
    trace
end

function tr(A::Matrix{<:Matrix}, B::Matrix{<:Matrix})
    trace = zero(eltype(A[1, 1]))
    @inbounds for j = 1:size(A, 1), i = 1:size(A, 2)
        trace += tr(A[i, j], B[j, i])
    end
    trace
end

