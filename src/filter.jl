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

    P = BlockDiagonalMatrix(joinpath(path, "covariance-point-sources-observed"), mmax,
                            progressbar=true, distribute=true)
    G = BlockDiagonalMatrix(joinpath(path, "covariance-galactic-observed"), mmax,
                            progressbar=true, distribute=true)
    S = BlockDiagonalMatrix(joinpath(path, "covariance-fiducial-signal-observed"), mmax,
                            progressbar=true, distribute=true)

    @. P = transfermatrix * point_sources * T(transfermatrix)
    @. G = transfermatrix * galactic      * T(transfermatrix)
    @. S = transfermatrix * signal        * T(transfermatrix)

    F = BlockDiagonalMatrix(joinpath(path, "complete-foregrounds-observed"), mmax,
                            progressbar=true, distribute=true)
    @. F = P + G

    # compute eig(F, S)
    # (this seems to be more numerically stable when F is the first argument)
    # apply the filter

    # N = noisematrix + F
    # compute eig(S, N)
    # (not sure which one should go first here...)
end



function whiten(transfermatrix, signal, noise)
end



