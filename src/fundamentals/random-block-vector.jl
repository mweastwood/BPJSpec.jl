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

struct WhiteNoiseBlockVector end
@inline Base.getindex(vector::WhiteNoiseBlockVector, idx...) = WhiteNoise()

struct RandomBlockVector{B}
    covariance :: B
end

function Base.getindex(vector::RandomBlockVector, idx...)
    # TODO: cache the results of the Cholesky decomposition
    C = vector.covariance[idx...]
    N = size(C, 1)
    U = chol(Hermitian(C))
    x = complex.(randn(N), randn(N)) ./ √2
    U'*x # random variable with the correct covariance
end

