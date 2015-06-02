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

function planewave(u,v,w;lmax::Int=100,mmax::Int=100)
    b = sqrt(u*u+v*v+w*w)
    θ = acos(w/b)
    ϕ = atan2(v,u)
    realpart = Alm(Complex128,lmax,mmax)
    imagpart = Alm(Complex128,lmax,mmax)
    for m = 0:mmax, l = m:lmax
        alm1 = 4π*(1im)^l*j(l,2π*b)*conj(Y(l,+m,θ,ϕ))
        alm2 = 4π*(1im)^l*j(l,2π*b)*conj(Y(l,-m,θ,ϕ))
        realpart[l,m] = (alm1 + conj(alm2))/2
        imagpart[l,m] = (alm1 - conj(alm2))/2im
    end
    realpart,imagpart
end

