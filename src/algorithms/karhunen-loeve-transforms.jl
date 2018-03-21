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
    foreground_filter!(output_mmodes, output_transfermatrix, output_covariance,
                       output_foreground_filter, output_noise_whitener,
                       input_mmodes, input_transfermatrix, input_noisematrix,
                       input_signalmatrix, input_foregroundmatrix;
                       threshold=0.1, tempdir="/tmp", cleanup=true)

Filter foreground emission from the dataset *and* whiten the noise.

**Arguments:**

* `output_mmodes` will be populated with the foreground filtered $m$-modes
* `output_transfermatrix` will be populated with the foreground-filtered transfer-matrix
* `output_covariance` will be populated with the covariance matrix of the output $m$-modes
* `output_foreground_filter` is the matrix used to filter foreground emission
* `output_noise_whitener` is the matrix used to whiten the noise covariance (the thermal noise and
  the foreground contribution)
* `input_mmodes` is the measured input $m$-modes
* `input_transfermatrix` is the transfer matrix describing the response of the interferometer
* `input_noisematrix` is the thermal noise covariance matrix
* `input_signalmatrix` is the covariance matrix due to the 21-cm power spectrum
* `input_foregroundmatrix` is the covariance matrix due to the foreground radio emission

!!! warn
    Please pay attention to the argument order. Swapping arguments could lead to this function
    filtering all of the 21-cm signal while leaving the foreground emission in tact!

**Keyword Arguments:**

* `threshold` is the maximum allowed foreground-signal ratio
* `tempdir` is a path to a directory where temporary files can be placed
* `cleanup` determines whether or not files placed in the `tempdir` are automatically removed. This
  can be set to `false` if you wish to check some of the intermediate output.
"""
function foreground_filter!(output_mmodes, output_transfermatrix, output_covariance,
                            output_foreground_filter, output_noise_whitener,
                            input_mmodes,  input_transfermatrix,  input_noisematrix,
                            input_signalmatrix, input_foregroundmatrix;
                            threshold=0.1, tempdir="/tmp", cleanup=true)
    mmax = input_mmodes.mmax
    storage(name) = MultipleFiles(joinpath(tempdir, name))

    # temporary matrices
    F  = create(MBlockMatrix, storage("observed-foreground-matrix"), mmax) |> ProgressBar
    S  = create(MBlockMatrix, storage("observed-signal-matrix"),     mmax) |> ProgressBar
    S′ = create(MBlockMatrix, storage("filtered-signal-matrix"),     mmax) |> ProgressBar
    N′ = create(MBlockMatrix, storage("filtered-noise-matrix"),      mmax) |> ProgressBar

    double_kl_transform!(ProgressBar(output_mmodes),
                         ProgressBar(output_transfermatrix),
                         ProgressBar(output_covariance),
                         ProgressBar(output_foreground_filter),
                         ProgressBar(output_noise_whitener),
                         input_mmodes, input_transfermatrix, input_noisematrix,
                         input_signalmatrix, input_foregroundmatrix,
                         F, S, S′, N′, threshold)

    if cleanup
        rm_old_blocks!(unwrap(F).storage)
        rm_old_blocks!(unwrap(S).storage)
        rm_old_blocks!(unwrap(S′).storage)
        rm_old_blocks!(unwrap(N′).storage)
    end

    output = (output_mmodes, output_transfermatrix, output_covariance,
              output_foreground_filter, output_noise_whitener)
    output
end

function double_kl_transform!(output_mmodes, output_transfermatrix, output_covariance,
                              output_foreground_filter, output_noise_whitener,
                              input_mmodes,  input_transfermatrix,  input_noisematrix,
                              input_signalmatrix, input_foregroundmatrix,
                              F, S, S′, N′, threshold)

    # Run the sky covariances through the interferometer's response.
    @. F = fix(input_transfermatrix * input_foregroundmatrix * T(input_transfermatrix))
    @. S = fix(input_transfermatrix *     input_signalmatrix * T(input_transfermatrix))
    N = input_noisematrix

    # Compute the foreground filter...
    V = output_foreground_filter
    construct_filter′(F, S) = construct_filter(F, S, threshold)
    @. V = construct_filter′(F, S)

    # ...and apply that filter.
    @. S′ = fix(T(V) * S * V)
    @. N′ = fix(T(V) * (N + F) * V)

    # Whiten the noise.
    W = output_noise_whitener
    @. W = construct_whiten(S′, N′)

    # Compute the final m-modes, transfer matrix, and covariance matrix.
    @. output_mmodes         = T(W) * (T(V) * input_mmodes)
    @. output_transfermatrix = T(W) * (T(V) * input_transfermatrix)
    @. output_covariance     = fix(T(W) * (S′ + N′) * W)

    nothing
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

