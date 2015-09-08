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

"""
    block(alm::Alm,m)

Get all the coefficients corresponding to a single value of m.
"""
function block(alm::Alm,m)
    a = zeros(Complex128,lmax(alm)-m+1)
    idx = 1
    for l = m:lmax(alm)
        a[idx] = alm[l,m]
        idx += 1
    end
    a
end

"""
    tikhonov(B::TransferMatrix,v::MModes;tolerance=0.0) -> Alm

Solve for the spherical harmonic coefficients with the given
transfer matrix and m-modes. `tolerance` acts as Tikhonov
regularization parameter that alleviates some of the numerical
problems with this process.
"""
function tikhonov(B::TransferMatrix,v::MModes;tolerance::Float64=0.0)
    alm = Alm(Complex128,lmax(B),mmax(B))
    for m = 0:mmax(B)
        Bm = B[m].block
        vm = v[m].block
        am = (Bm'*Bm + tolerance*I)\Bm'*vm
        idx = 1
        for l = m:lmax(B)
            alm[l,m] = am[idx]
            idx += 1
        end
    end
    alm
end

