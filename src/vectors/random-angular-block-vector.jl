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

struct RandomAngularBlockVector <: BlockVector
    lmax :: Int
    mmax :: Int
    covariance :: AngularCovarianceMatrix
end

function RandomAngularBlockVector(covariance::AngularCovarianceMatrix)
    lmax = mmax = covariance.lmax
    RandomAngularBlockVector(lmax, mmax, covariance)
end

indices(matrix::RandomAngularBlockVector) =
    [(l, m) for m = 0:matrix.mmax for l = m:matrix.lmax]

function Base.getindex(vector::RandomAngularBlockVector, l, m)
    C = vector.covariance[l, m]
    N = size(C, 1)
    U = chol(C)
    x = complex.(randn(N), randn(N)) ./ √2
    U'*x # random variable with the correct covariance
end

