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

# Design decision
# ===============
# We are going to be compressing the transfer matrix using its singular value decomposition, however
# we'd like to also use the results of the SVD to compress the m-modes. The essential problem
# however is that the SVD is just as large as the transfer matrix itself, so it will take a lot of
# disk space to store. Therefore we will compute the SVD, compress both the transfer matrix and the
# m-modes before discarding the SVD.

function full_rank_compress(transfermatrix, noisematrix)
    path = dirname(transfermatrix.path)
    mmax = transfermatrix.mmax
    frequencies = transfermatrix.frequencies
    bandwidth   = transfermatrix.bandwidth

    suffix = "-compressed"
    file = joinpath(path, "transfer-matrix"*suffix)
    output_transfermatrix = SpectralBlockDiagonalMatrix(file, mmax, frequencies, bandwidth,
                                                        progressbar=true, distribute=true)
    file = joinpath(path, "noise-matrix"*suffix)
    output_noisematrix = SpectralBlockDiagonalMatrix(file, mmax, frequencies, bandwidth,
                                                     progressbar=true, distribute=true)

    multi_broadcast!(_full_rank_compress, compress,
                     (output_transfermatrix, output_noisematrix),
                     (transfermatrix, noisematrix))
end

function _full_rank_compress(B, N)
    F = svdfact(B)
    B′ = F[:U]'*B
    N′ = F[:U]'*N*F[:U]
    B′, N′
end

