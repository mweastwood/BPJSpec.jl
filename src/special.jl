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

import GSL: sf_legendre_sphPlm

"""
The spherical harmonic function.
"""
function Y(l,m,θ,ϕ)
    sf_legendre_sphPlm(l,abs(m),cos(θ))*exp(1im*m*ϕ)
end

"""
The spherical Bessel function (of the first kind).
"""
function j(l,x)
    sqrt(π/(2x))*besselj(l+1/2,x)
end

"""
The trace of the product of two matrices, computed
without first computing the product itself.

    tr(AB)
"""
function tr{T}(A::Matrix{T},B::Matrix{T})
    trace = zero(T)
    Bt = transpose(B)
    for i in eachindex(A)
        trace += A[i]*Bt[i]
    end
    trace
end

