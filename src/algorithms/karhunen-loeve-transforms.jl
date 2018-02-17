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

function kltransforms(transfermatrix, noisematrix, foregroundmatrix, signalmatrix)
    mmax = transfermatrix.mmax

    # Run the sky covariances through the interferometer's response.
    file = joinpath(dirname(foregroundmatrix.path), "covariance-matrix-fiducial-foregrounds-observed")
    F = DenseBlockDiagonalMatrix(file, mmax, progressbar=true, distribute=true)
    file = joinpath(dirname(signalmatrix.path), "covariance-matrix-fiducial-signal-observed")
    S = DenseBlockDiagonalMatrix(file, mmax, progressbar=true, distribute=true)
    @. F = fix(transfermatrix * foregroundmatrix * T(transfermatrix))
    @. S = fix(transfermatrix *     signalmatrix * T(transfermatrix))

    # Compute the foreground filter...
    file = joinpath(dirname(transfermatrix.path), "foreground-filter")
    V = DenseBlockDiagonalMatrix(file, mmax, progressbar=true, distribute=true)
    @. V = construct_filter(F, S)

    # ...and apply that filter.
    file = joinpath(dirname(signalmatrix.path), "covariance-matrix-fiducial-signal-filtered")
    S′ = DenseBlockDiagonalMatrix(file, mmax, progressbar=true, distribute=true)
    file = joinpath(dirname(noisematrix.path), "covariance-matrix-noise-filtered")
    N′ = DenseBlockDiagonalMatrix(file, mmax, progressbar=true, distribute=true)
    @. S′ = fix(T(V) * S * V)
    @. N′ = fix(T(V) * noisematrix * V) + H(T(V) * F * V)

    # Whiten the noise.
    file = joinpath(dirname(transfermatrix.path), "whitening-matrix")
    W = DenseBlockDiagonalMatrix(file, mmax, progressbar=true, distribute=true)
    @. W = construct_whiten(S′, N′)

    # Compute the final transfer matrix and covariance matrix.
    file = joinpath(dirname(transfermatrix.path), "transfer-matrix-final")
    B = DenseBlockDiagonalMatrix(file, mmax, progressbar=true, distribute=true)
    file = joinpath(dirname(transfermatrix.path), "covariance-matrix-final")
    C = DenseBlockDiagonalMatrix(file, mmax, progressbar=true, distribute=true)
    @. B = T(W) * (T(V) * transfermatrix)
    @. C = H(T(W) * S′ * W) + H(T(W) * N′ * W)

    B, C
end

function construct_filter(F, S)
    # We want to make sure that we're hitting the correct branch in `eig` for Hermitian, positive
    # definite matrices. So let's verify that here and make a noisy error if it's not true.
    isposdef(F) || error("F is not positive definite")
    isposdef(S) || error("S is not positive definite")
    λ, V = eig(F, S) # Note: this seems to be more numerically stable than eig(S, F)
    idx = searchsortedlast(λ, 0.1) # foreground-signal ratio < 0.1
    V[:, 1:idx]
end

function construct_whiten(S, N)
    isposdef(S) || error("S is not positive definite")
    isposdef(N) || error("N is not positive definite")
    λ, V = eig(S, N)
    V
end

