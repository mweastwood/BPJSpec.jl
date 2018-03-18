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
    full_rank_compress!(output_mmodes, output_transfermatrix, output_noisematrix,
                        input_mmodes,  input_transfermatrix,  input_noisematrix;
                        progress=false)

In the case where the interferometer has more baselines than there are spherical harmonic
coefficients to measure, the transfer matrix is tall and skinny. This also indicates that we have
made redundant measurements that can be averaged together with no information loss.

In this routine we use the singular value decomposition (SVD) of the transfer matrix to compress the
measurements. However, the SVD is just as large as the transfer matrix itself, and will take a lot
of disk space to store. Therefore we will compute the SVD, compress everything with it all at once
so that there is no need to store the SVD as well.

**Arguments:**

* `output_mmodes` the output compressed $m$-modes
* `output_transfermatrix` the output compressed transfer matrix
* `output_noisematrix` the output compressed noise covariance matrix
* `input_mmodes` the input $m$-modes that will be compressed
* `input_transfermatrix` the input transfer matrix that will be used to generate the compression
* `input_noisematrix` the input noise covariance matrix

**Keyword Arguments:**

* `progress` if set to `true`, a progress bar will be displayed
"""
function full_rank_compress!(output_mmodes, output_transfermatrix, output_noisematrix,
                             input_mmodes,  input_transfermatrix,  input_noisematrix;
                             progress=false)
    multi_broadcast!(_full_rank_compress,
                     (output_mmodes, output_transfermatrix, output_noisematrix),
                     (input_mmodes,  input_transfermatrix,  input_noisematrix),
                     progress=progress)
end

function _full_rank_compress(v, B, N)
    F = svdfact(B)
    U = F[:U]
    v′ = U'*v
    B′ = U'*B
    N′ = fix(U'*N*U)
    v′, B′, N′
end

