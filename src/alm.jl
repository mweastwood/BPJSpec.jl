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

function block(alm::Alm,m)
    a = zeros(Complex128,lmax(alm)-m+1)
    idx = 1
    for l = m:lmax(alm)
        a[idx] = alm[l,m]
        idx += 1
    end
    a
end

function Alm(transfermatrix::TransferMatrix,mmodes::MModes;
             regularization_parameter::Float64=1.0)
    alm = Alm(Complex128,lmax(transfermatrix),
                         mmax(transfermatrix))
    for m = 0:mmax(transfermatrix)
        B = block(transfermatrix,m)
        v = block(mmodes,m)
        a = (B'*B + regularization_parameter*I)\B'*v
        idx = 1
        for l = m:lmax(transfermatrix)
            alm[l,m] = a[idx]
            idx += 1
        end
    end
    alm
end

\(transfermatrix::TransferMatrix,mmodes::MModes) = Alm(transfermatrix,mmodes)

