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

function noise_bias(transfermatrix, covariancematrix, basis; iterations=10)
    N = length(basis)
    M = length(workers())
    b = zeros(N)
    @sync for worker in workers()
        @async begin
            b′ = remotecall_fetch(_bias, worker, transfermatrix, covariancematrix,
                                  basis, cld(iterations, M))
            b .+= b′
        end
    end
    b ./= M
    b
end

function _bias(transfermatrix, covariancematrix, basis, iterations)
    cache!(transfermatrix)
    cache!(covariancematrix)
    foreach(cache!, basis)
    compute_bias(transfermatrix, covariancematrix, basis, iterations)
end

function compute_bias(B, C, basis, iterations)
    N = length(basis)
    lmax = mmax = B.mmax

    w   = WhiteNoiseVector()
    Bw  = create(MBlockVector, mmax)
    CBw = create(AngularBlockVector, lmax, mmax)

    b = zeros(N)
    q = zeros(N)

    for iter = 1:iterations
        fisher_iteration!(q, B, C, basis, w, Bw, CBw)
        b .+= q
    end
    b ./= iterations
    b
end

