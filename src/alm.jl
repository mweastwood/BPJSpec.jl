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

"""
    tikhonov(B::TransferMatrix, v::MModes, tolerance)

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
function tikhonov(B::TransferMatrix, v::MModes, tolerance)
    # note for now we assume that B and v each only have one frequency channel
    alm = Alm(Complex128, B.lmax, B.mmax)
    p = Progress(B.mmax+1, "Solving: ")
    for m = 0:B.mmax
        Bm = B[m,1]
        vm = v[m,1]
        account_for_flags!(Bm, vm)
        am = tikhonov(Bm, vm, tolerance)
        setblock!(alm, am, m)
        next!(p)
    end
    alm
end

function tikhonov(A::Matrix, b::Vector, tol)
    At = A'
    AA = At*A # rate limiting step
    Ab = At*b
    D = tol*I
    out = (AA + D)\Ab
    out
end

"""
    account_for_flags!(transfermatrix_block, mmodes_block)

Flagged baselines will have their m-mode equal to zero. We therefore
need to set the corresponding row in the transfer matrix to zero as well.
"""
function account_for_flags!(transfermatrix_block, mmodes_block)
    for α = 1:length(mmodes_block)
        if abs(mmodes_block[α]) < eps(Float64)
            transfermatrix_block[α,:] = 0
        end
    end
end

