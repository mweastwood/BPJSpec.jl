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

function full_rank_compress(input_mmodes, input_transfermatrix, input_noisematrix;
                            mmodes_storage=NoFile(),
                            transfermatrix_storage=NoFile(),
                            noisematrix_storage=NoFile(),
                            progress=false)

    output_mmodes         = similar(input_mmodes,                 mmodes_storage)
    output_transfermatrix = similar(input_transfermatrix, transfermatrix_storage)
    output_noisematrix    = create(MFBlockMatrix, noisematrix_storage, input_noisematrix.mmax,
                                   input_noisematrix.frequencies, input_noisematrix.bandwidth)

    multi_broadcast!(_full_rank_compress,
                     (output_mmodes, output_transfermatrix, output_noisematrix),
                     (input_mmodes, input_transfermatrix, input_noisematrix), progress=progress)
end

function _full_rank_compress(v, B, N)
    F = svdfact(B)
    v′ = F[:U]'*v
    B′ = F[:U]'*B
    N′ = fix(F[:U]'*N*F[:U])
    v′, B′, N′
end

