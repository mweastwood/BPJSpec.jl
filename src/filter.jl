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

function filter(transfermatrix, noisematrix)
    path = dirname(transfermatrix.path)
    lmax = mmax = transfermatrix.mmax
    frequencies = transfermatrix.frequencies

    point_sources = AngularCovarianceMatrix(joinpath(path, "covariance-point-sources"),
                                            lmax, frequencies, extragalactic_point_sources())
    galactic      = AngularCovarianceMatrix(joinpath(path, "covariance-galactic"),
                                            lmax, frequencies, galactic_synchrotron())
    signal        = AngularCovarianceMatrix(joinpath(path, "covariance-fiducial-signal"),
                                            lmax, frequencies, fiducial_signal_model())

    # Run the sky covariances through the interferometer's response.
    P = BlockDiagonalMatrix(joinpath(path, "covariance-point-sources-observed"), mmax,
                            progressbar=true, distribute=true)
    G = BlockDiagonalMatrix(joinpath(path, "covariance-galactic-observed"), mmax,
                            progressbar=true, distribute=true)
    S = BlockDiagonalMatrix(joinpath(path, "covariance-fiducial-signal-observed"), mmax,
                            progressbar=true, distribute=true)
    F = BlockDiagonalMatrix(joinpath(path, "covariance-complete-foregrounds-observed"), mmax,
                            progressbar=true, distribute=true)
    @. P = H(transfermatrix * point_sources * T(transfermatrix))
    @. G = H(transfermatrix * galactic      * T(transfermatrix))
    @. S = H(transfermatrix * signal        * T(transfermatrix))
    @. F = P + G

    # Compute the foreground filter...
    V = BlockDiagonalMatrix(joinpath(path, "foreground-filter"), mmax,
                            progressbar=true, distribute=true)
    @. V = construct_filter(F, S)

    # ...and apply that filter.
    S′ = BlockDiagonalMatrix(joinpath(path, "covariance-fiducial-signal-filtered"), mmax,
                             progressbar=true, distribute=true)
    N′ = BlockDiagonalMatrix(joinpath(path, "covariance-compelete-noise-filtered"), mmax,
                             progressbar=true, distribute=true)
    @. S′ = H(T(V) * S * V)
    @. N′ = H(T(V) * N * V) + H(T(V) * F * V)

    # Whiten the noise.
    # ...
end

function construct_filter(F, S)
    λ, V = eig(F, S) # Note: this seems to be more numerically stable than eig(S, F)
    idx = searchsortedlast(λ, 0.1) # foreground-signal ratio < 0.1
    V[:, 1:idx]
end

