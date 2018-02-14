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

function fisher(transfermatrix, covariancematrix, basis)
    @time cache!(transfermatrix)
    @time cache!(covariancematrix)
    @time foreach(cache!, basis)
    save("/dev/shm/mweastwood/cache.jld2",
         "B", transfermatrix, "C", covariancematrix, "Ca", basis)
    #@time transfermatrix, covariancematrix, basis =
    #    load("/dev/shm/mweastwood/cache.jld2", "B", "C", "Ca")
    @time _fisher(transfermatrix, covariancematrix, basis)
end

function _fisher(B, C, basis, iterations=10)
    N = length(basis)
    lmax = mmax = B.mmax

    v   = RandomVector(C)
    Bv  = BlockDiagonalVector(mmax)
    CBv = AngularBlockVector(lmax, mmax)

    prg = Progress(iterations)
    q = zeros(N, iterations)
    for iter in 1:iterations
        fisher_iteration!(view(q, :, iter), B, C, basis, v, Bv, CBv)
        next!(prg)
    end
    q
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

