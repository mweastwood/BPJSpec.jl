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






#struct RandomVector{B}
#    covariance :: B
#end
#
#function Base.getindex(vector::RandomVector, idx)
#    C = vector.covariance[idx]
#    N = size(C, 1)
#    U = chol(C)
#    x = complex.(randn(N), randn(N)) ./ √2
#    U'*x # random variable with the correct covariance
#end
#
#struct WhiteNoiseVector  end
#struct WhiteNoiseMachine end
#
#@inline Base.getindex(vector::WhiteNoiseVector, idx) = WhiteNoiseMachine()
#@inline (machine::WhiteNoiseMachine)() = complex(randn(), randn())/√2
#
#function Base.:*(A::Matrix, w::WhiteNoiseMachine)
#    N, M = size(A)
#    v = [w() for idx = 1:M]
#    A*v
#end
#
#function Base.:*(w::WhiteNoiseMachine, A::Matrix)
#    N, M = size(A)
#    v = [w() for idx = 1:N]
#    v'*A
#end
#
#function Base.:/(w::WhiteNoiseMachine, A::Matrix)
#    N, M = size(A)
#    v = [w() for idx = 1:M]
#    v'/A
#end
#
#function Base.:+(v::Vector, w::WhiteNoiseMachine)
#    [element+w() for element in v]
#end
#Base.:+(w::WhiteNoiseMachine, v::Vector) = v + w
#
#Base.transpose(w::WhiteNoiseMachine) = w

