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

function fisher_information(transfermatrix, covariancematrix, basis; iterations=10)
    N = length(basis)
    M = length(workers())
    F = zeros(N, N)
    @sync for worker in workers()
        @async begin
            F′ = remotecall_fetch(_fisher, worker, transfermatrix, covariancematrix,
                                  basis, cld(iterations, M))
            F .+= F′
        end
    end
    F ./= M
    F
end

function _fisher(transfermatrix, covariancematrix, basis, iterations)
    cache!(transfermatrix)
    cache!(covariancematrix)
    foreach(cache!, basis)
    compute_fisher(transfermatrix, covariancematrix, basis, iterations)
end

function compute_fisher(B, C, basis, iterations)
    N = length(basis)
    lmax = mmax = B.mmax

    v   = RandomVector(C)
    Bv  = BlockDiagonalVector(mmax)
    CBv = AngularBlockVector(lmax, mmax)

    q  = zeros(N)
    μ  = zeros(N)
    σ² = zeros(N, N)

    for iter = 1:iterations
        fisher_iteration!(q, B, C, basis, v, Bv, CBv)
        μ  .+= q
        σ² .+= q*q'
    end
    μ  ./= iterations
    σ² ./= iterations
    σ² .- μ*μ'
end

function fisher_iteration!(q, B, C, basis, v, Bv, CBv)
    N = length(basis)
    @. Bv = T(B) * (C \ v)
    Bv′ = AngularBlockVector(Bv)
    for a = 1:N
        Ca = basis[a]
        @. CBv = Ca * Bv′
        q[a] = real(dot(CBv, Bv′))
    end
    q
end

