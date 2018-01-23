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

"""
    Cforeground(l, ν1, ν2, ν0, A, α, β, ζ)

Evaluate a model for the multifrequency angular power spectrum of a foreground component. Note that
`ν0`, `A`, `α`, `β`, and `ζ` are all parameters of the model.
"""
function Cforeground(l, ν1, ν2, ν0, A, α, β, ζ)
    (A * (l+1)^(-α)
        * (ν1*ν2/ν0^2)^(-β)
        * exp(-log(ν1/ν2)^2/(2*ζ^2)))
end

immutable ForegroundComponent <: SkyComponent
    ν0 :: Float64
    A  :: Float64
    α  :: Float64
    β  :: Float64
    ζ  :: Float64
end

function (::ForegroundComponent)(l, ν1, ν2)
    Cforeground(l, ν1, ν2, model.ν0, model.A, model.α, model.β, model.ζ)
end

