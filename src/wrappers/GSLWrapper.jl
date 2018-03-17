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

module GSLWrapper

export Y

using GSL

"""
    Y(l, m, θ, ϕ)

The spherical harmonic function (using the Condon-Shortley phase convention).

**Usage:**

```jldoctest
julia> using BPJSpec

julia> BPJSpec.Y(0, 0, 1, 2)
0.28209479177387814 + 0.0im

julia> 1/sqrt(4π)
0.28209479177387814
```
"""
function Y(l, m, θ, ϕ)
    out = GSL.sf_legendre_sphPlm(l, abs(m), cos(θ)) * cis(m*ϕ)
    # GSL already applies the Condon-Shortley phase convention to this result, so if we want m < 0
    # we need to undo the factor of (-1)^m.
    out *= m < 0? (-1)^(-m) : 1
    out
end

end

