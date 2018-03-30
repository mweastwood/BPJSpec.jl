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
    tikhonov(transfermatrix, mmodes; regularization=1e-2, mfs=false, storage=NoFile())

Create a dirty image of the sky using Tikhonov regularization.

**Arguments:**

* `transfermatrix` the interferometer's transfer matrix, describing its response to the sky
* `mmodes` the $m$-modes measured by the interferometer

**Keyword Arguments:**

* `regularization` the amplitude of the Tikhonov regularization parameter
* `mfs` determines whether or not to perform Multi-Frequency Synthesis imaging. If this parameter is
  set to `true`, all frequency channels will be used to generate a single image of the sky. If this
  parameter is set to `false`, an image of the sky will be generated for each frequency channel.
* `storage` determines how the computed spherical harmonic coefficients will be stored
"""
function tikhonov(transfermatrix, mmodes; regularization=1.0, mfs=false, storage=NoFile())
    if mfs
        return tikhonov_mfs(transfermatrix, mmodes, regularization, storage)
    else
        return tikhonov_nothing_special(transfermatrix, mmodes, regularization, storage)
    end
end

function tikhonov_nothing_special(transfermatrix, mmodes, regularization, storage)
    alm = similar(mmodes, storage)
    invert(B, v) = _tikhonov_nothing_special(B, v, regularization)
    @. alm = invert(transfermatrix, mmodes)
    alm
end

function _tikhonov_nothing_special(B, v, regularization)
    BLAS.set_num_threads(16)
    B, v = _tikhonov_propagate_flags(B, v)
    _tikhonov_inversion(B'*B, B'*v, regularization)
end

function _tikhonov_propagate_flags(B, v)
    f = v .== 0 # flags
    v = v[.!f]
    B = B[.!f, :]
    B, v
end

function _tikhonov_inversion(BB, Bv, regularization)
    (BB + regularization*I) \ Bv
end

function tikhonov_mfs(transfermatrix, mmodes, regularization, storage)
    lmax = mmax = transfermatrix.mmax
    alm  = create(MBlockVector, storage, mmax)

    # Try to make the scaling of the regularization parameter consistent between regular and
    # multi-frequency synthesis images
    regularization *= length(transfermatrix.frequencies)

    pool  = CachingPool(workers())
    queue = collect(0:mmax)

    lck = ReentrantLock()
    prg = Progress(length(queue))
    increment() = (lock(lck); next!(prg); unlock(lck))

    @sync for worker in workers()
        @async while length(queue) > 0
            m = shift!(queue)
            alm[m] = remotecall_fetch(_tikhonov_mfs, pool, transfermatrix, mmodes,
                                      regularization, lmax, m)
            increment()
        end
    end
    alm
end

function _tikhonov_mfs(transfermatrix, mmodes, regularization, lmax, m)
    BLAS.set_num_threads(16)
    BB = zeros(Complex128, lmax-m+1, lmax-m+1)
    Bv = zeros(Complex128, lmax-m+1)
    for β = 1:length(mmodes.frequencies)
        _tikhonov_accumulate!(BB, Bv, transfermatrix[m, β], mmodes[m, β])
    end
    _tikhonov_inversion(BB, Bv, regularization)
end

function _tikhonov_accumulate!(BB, Bv, B, v)
    B, v = _tikhonov_propagate_flags(B, v)
    BB .+= B'*B
    Bv .+= B'*v
end

