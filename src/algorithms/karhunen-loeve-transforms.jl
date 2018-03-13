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

function kltransforms(mmodes, transfermatrix, noisematrix,
                      foregroundmatrix, signalmatrix,
                      observed_foregroundmatrix_storage = NoFile(),
                      observed_signalmatrix_storage     = NoFile(),
                      foreground_filter_storage         = NoFile(),
                      filtered_signalmatrix_storage     = NoFile(),
                      filtered_noisematrix_storage      = NoFile(),
                      whitening_matrix_storage          = NoFile(),
                      mmodes_storage                    = NoFile(),
                      transfermatrix_storage            = NoFile(),
                      covariancematrix_storage          = NoFile();
                      threshold = 0.1)
    mmax = mmodes.mmax

    # Run the sky covariances through the interferometer's response.
    F = create(MBlockMatrix, observed_foregroundmatrix_storage, mmax) |> ProgressBar
    S = create(MBlockMatrix, observed_signalmatrix_storage,     mmax) |> ProgressBar
    @. F = fix(transfermatrix * foregroundmatrix * T(transfermatrix))
    @. S = fix(transfermatrix *     signalmatrix * T(transfermatrix))

    # Compute the foreground filter...
    V = create(MBlockMatrix, foreground_filter_storage, mmax) |> ProgressBar
    construct_filter′(F, S) = construct_filter(F, S, threshold)
    @. V = construct_filter′(F, S)

    # ...and apply that filter.
    S′ = create(MBlockMatrix, filtered_signalmatrix_storage, mmax) |> ProgressBar
    N′ = create(MBlockMatrix, filtered_noisematrix_storage,  mmax) |> ProgressBar
    @. S′ = fix(T(V) * S * V)
    @. N′ = fix(T(V) * noisematrix * V) + H(T(V) * F * V)

    # Whiten the noise.
    W = create(MBlockMatrix, whitening_matrix_storage, mmax) |> ProgressBar
    @. W = construct_whiten(S′, N′)

    # Compute the final m-modes, transfer matrix, and covariance matrix.
    v = create(MBlockVector, mmodes_storage,           mmax) |> ProgressBar
    B = create(MBlockMatrix, transfermatrix_storage,   mmax) |> ProgressBar
    C = create(MBlockMatrix, covariancematrix_storage, mmax) |> ProgressBar
    @. v = T(W) * (T(V) * mmodes)
    @. B = T(W) * (T(V) * transfermatrix)
    @. C = H(T(W) * S′ * W) + H(T(W) * N′ * W)

    unwrap(v), unwrap(B), unwrap(C)
end

function construct_filter(F, S, threshold)
    # We want to make sure that we're hitting the correct branch in `eig` for Hermitian, positive
    # definite matrices. So let's verify that here and make a noisy error if it's not true.
    isposdef(F) || error("F is not positive definite")
    isposdef(S) || error("S is not positive definite")
    λ, V = eig(F, S) # Note: this seems to be more numerically stable than eig(S, F)
    idx = searchsortedlast(λ, threshold) # foreground-signal ratio < threshold
    V[:, 1:idx]
end

function construct_whiten(S, N)
    isposdef(S) || error("S is not positive definite")
    isposdef(N) || error("N is not positive definite")
    λ, V = eig(S, N)
    V
end

