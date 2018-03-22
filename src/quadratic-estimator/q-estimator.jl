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

doc"""
    q_estimator(mmodes, transfermatrix, covariancematrix, basis)

Evaluate the $q$ estimator:

```math
q_a = v^* C^{-1} B C_a B^* C^{-1} v
```

**Arguments:**

* `mmodes` or $v$ specifies the list of measured $m$-modes
* `transfermatrix` or $B$ specifies the interferometer's response to the sky
* `covariancematrix` or $C$ specifies the covariance of the measured $m$-modes
* `basis` or $C_a$ is a list of angular covariance matrices that represent the change in the
  covariance with respect to an increase in power of each 21-cm power spectrum bin
"""
function q_estimator(mmodes, transfermatrix, covariancematrix, basis)
    N = length(basis)
    output = zeros(N)
    q_estimator!(output, mmodes, transfermatrix, covariancematrix, basis)
    output
end

function q_estimator!(output, mmodes, transfermatrix, covariancematrix, basis)
    N = length(basis)
    lmax = mmax = transfermatrix.mmax
    frequencies = first(basis).frequencies
    bandwidth   = first(basis).bandwidth

    Bv = create(MBlockVector, mmax)
    @. Bv = T(transfermatrix) * (covariancematrix \ mmodes)
    Bv′ = create(MultiFrequencyAlm, Bv, frequencies, bandwidth)
    CBv = create(LMBlockVector, lmax, mmax, frequencies, bandwidth)
    for a = 1:N
        Ca = basis[a]
        @. CBv = Ca * Bv′
        output[a] = real(dot(CBv, Bv′))
    end
    output
end

