# Copyright (c) 2015 Michael Eastwood
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

# For the most part HEALPix.Alm will be sufficient for our purposes,
# however we will need some additional definitions to make things
# consistent with the language of m-mode analysis.

doc"""
    getblock(alm::Alm, m)

Get all the coefficients corresponding to a single value of $m$.
"""
function getblock(alm::Alm, m)
    a = zeros(Complex128,lmax(alm)-m+1)
    idx = 1
    for l = m:lmax(alm)
        a[idx] = alm[l,m]
        idx += 1
    end
    a
end

doc"""
    setblock!(alm::Alm, x, m)

Set all the coefficients corresponding to a single value of $m$
to a value given by `x`.
"""
function setblock!(alm::Alm, x, m)
    idx = 1
    for l = m:lmax(alm)
        alm[l,m] = x[idx]
        idx += 1
    end
    alm
end

function *(B::TransferMatrix,alm::Alm)
    is_single_frequency(B) || error("Expected single-frequency transfer matrix.")
    lmax(B) == lmax(alm) || error("The values of lmax must be the same.")
    mmax(B) == mmax(alm) || error("The values of mmax must be the same.")
    blocks = VectorBlock[]
    for m = 0:mmax(B)
        Bm = B[m+1].block
        am = getblock(alm,m)
        vm = Bm*am
        push!(blocks,VectorBlock(vm))
    end
    meta = MModesMeta(B.meta.m,B.meta.ν)
    MModes(blocks,meta)
end

"""
    tikhonov(B::TransferMatrix, v::MModes; tolerance=0.0) -> Alm

Solve for the spherical harmonic coefficients with the given
transfer matrix and m-modes. `tolerance` acts as Tikhonov
regularization parameter that alleviates some of the numerical
problems with this process.

Note that the interferometer cannot see the entire sky. A northern hemisphere
telescope will never see the southern hemisphere. This means that there
must be some information about the sky that is lost during the act of
observing the sky. In linear algebra this means that the transfer matrix
must be singular, and we will therefore need some amount of numerical
regularization to invert it.
"""
function tikhonov(B::TransferMatrix, v::MModes;
                  tolerance::Float64 = 0.0)
    is_single_frequency(B) || error("Expected single-frequency transfer matrix.")
    is_single_frequency(v) || error("Expected single-frequency m-modes.")
    alm = Alm(Complex128,lmax(B),mmax(B))
    for m = 0:mmax(B)
        Bm = B[m+1].block
        vm = v[m+1].block
        am = (Bm'*Bm + tolerance*I)\Bm'*vm
        setblock!(alm,am,m)
    end
    alm
end

