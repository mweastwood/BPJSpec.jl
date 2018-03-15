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

struct WhiteNoise end
@inline (noise::WhiteNoise)(N=1) = complex.(randn(N), randn(N))/√2
@inline Base.getindex(noise::WhiteNoise, idx) = noise(length(idx))

Base.:*(A::AbstractMatrix, w::WhiteNoise) = ((N, M) = size(A); A*w(M))
Base.:\(A::AbstractMatrix, w::WhiteNoise) = ((N, M) = size(A); A\w(N))

Base.:*(v::AbstractVector, w::WhiteNoise) = (N = length(v); v*w(M))
Base.:+(v::AbstractVector, w::WhiteNoise) = (N = length(v); v+w(N))
Base.:+(w::WhiteNoise, v::AbstractVector) = v+w
Base.:-(v::AbstractVector, w::WhiteNoise) = (N = length(v); v-w(N))
Base.:-(w::WhiteNoise, v::AbstractVector) = -(v-w)

