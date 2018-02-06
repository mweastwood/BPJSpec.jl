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

function full_rank_compress(transfermatrix::TransferMatrix, noise::NoiseMatrix)
    path = dirname(transfermatrix.path)
    lmax = getlmax(input)
    mmax = lmax
    frequencies = transfermatrix.metadata.frequencies

    output_transfermatrix = SpectralBlockDiagonalMatrix(joinpath(path, "transfermatrix-compressed"),
                                                        mmax, transfermatrix.metadata.frequencies,
                                                        progressbar=true, distribute=true)

    output_noise = SpectralBlockDiagonalMatrix(joinpath(path, "noisematrix-compressed"),
                                               mmax, transfermatrix.metadata.frequencies,
                                               progressbar=true, distribute=true)

    multi_broadcast!(compress, (output_transfermatrix, output_noise), (transfermatrix, noise))
end

function compress(B, N)
    F = svdfact(B)
    B′ = F[:U]'*B
    N′ = F[:U]'*N*F[:U]
    B′, N′
end

function average_channels(transfermatrix::SpectralBlockDiagonalMatrix,
                          noise::SpectralBlockDiagonalMatrix, Navg)
    Nfreq  = length(input.frequencies)
    Nfreq′ = Nfreq ÷ Navg + 1
    mmax  = input.mmax

    frequencies = zeros(eltype(input.frequencies), Nfreq′)
    for idx = 1:Nfreq′
        range = (1:Navg) + (idx-1)*Navg
        range = range[1]:min(range[end], Nfreq)
        frequencies[idx] = mean(input.frequencies[range])
    end

    output = SpectralBlockDiagonalMatrix(path, mmax, frequencies)
    prg = Progress(Nfreq′)
    for idx = 1:Nfreq′
        range = (1:Navg) + (idx-1)*Navg
        range = range[1]:min(range[end], Nfreq)
        for m = 0:mmax
            block = input[m, range[1]]
            for β in range[2:end]
                block = [block; input[m, β]]
            end
            block = _lossless_compress_svd(block)
            output[m, frequencies[idx]] = block
        end
        next!(prg)
    end
end

