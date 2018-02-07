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

function full_rank_compress(transfermatrix::TransferMatrix, noisematrix::NoiseMatrix)
    path = dirname(transfermatrix.path)
    lmax = getlmax(transfermatrix)
    mmax = lmax
    frequencies = transfermatrix.metadata.frequencies

    suffix = "-compressed"
    file = joinpath(path, "transfer-matrix"*suffix)
    output_transfermatrix = SpectralBlockDiagonalMatrix(file, mmax, frequencies,
                                                        progressbar=true, distribute=true)
    file = joinpath(path, "noise-matrix"*suffix)
    output_noisematrix = SpectralBlockDiagonalMatrix(file, mmax, frequencies,
                                                     progressbar=true, distribute=true)

    multi_broadcast!(compress, (output_transfermatrix, output_noisematrix),
                     (transfermatrix, noisematrix))
end

function compress(B, N)
    F = svdfact(B)
    B′ = F[:U]'*B
    N′ = F[:U]'*N*F[:U]
    B′, N′
end

function average_channels(transfermatrix::SpectralBlockDiagonalMatrix,
                          noisematrix::SpectralBlockDiagonalMatrix, Navg)
    path = dirname(transfermatrix.path)
    Nfreq  = length(transfermatrix.frequencies)
    Nfreq′ = Nfreq ÷ Navg + 1
    mmax   = transfermatrix.mmax

    # Decide which channels to average down and the resulting average frequencies.
    ranges = Array{UnitRange{Int}}(Nfreq′)
    frequencies = Array{eltype(transfermatrix.frequencies)}(Nfreq′)
    for idx = 1:Nfreq′
        range = (1:Navg) + (idx-1)*Navg
        range = range[1]:min(range[end], Nfreq)
        ranges[idx] = range
        frequencies[idx] = mean(transfermatrix.frequencies[range])
    end

    suffix = "-averaged"
    file = joinpath(path, "transfer-matrix"*suffix)
    output_transfermatrix = SpectralBlockDiagonalMatrix(file, mmax, frequencies,
                                                        progressbar=true, distribute=false)
    file = joinpath(path, "noise-matrix"*suffix)
    output_noisematrix = SpectralBlockDiagonalMatrix(file, mmax, frequencies,
                                                     progressbar=true, distribute=false)

    # Perform the averaging.
    prg = Progress(Nfreq′)
    for idx = 1:Nfreq′
        range = ranges[idx]
        for m = 0:mmax
            B = stack(transfermatrix, m, range)
            N = noisematrix[m, range]
            B′, N′ = compress(B, N)
            output_transfermatrix[m, idx] = B′
            output_noisematrix[m, idx] = B′
        end
        next!(prg)
    end

    output_transfermatrix, output_noisematrix
end

