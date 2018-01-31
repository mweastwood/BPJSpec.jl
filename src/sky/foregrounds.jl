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
    x = (l+1)/1000
    y = uconvert(Unitful.NoUnits, ν1*ν2/ν0^2)
    z = log(ν1/ν2)
    ustrip(uconvert(u"K^2", A * x^-α * y^-β * exp(-z^2/(2ζ^2))))
end

struct ForegroundComponent <: SkyComponent
    ν0 :: typeof(1.0u"Hz")
    A  :: typeof(1.0u"K^2")
    α  :: Float64
    β  :: Float64
    ζ  :: Float64
end

function (component::ForegroundComponent)(l, ν1, ν2)
    Cforeground(l, ν1, ν2, component.ν0, component.A, component.α, component.β, component.ζ)
end

function extragalactic_point_sources()
    # http://adsabs.harvard.edu/abs/2005ApJ...625..575S
    ForegroundComponent(130.0u"MHz", 57.0u"mK^2", 1.1, 2.07, 1.0)
end

function extragalactic_free_free()
    # http://adsabs.harvard.edu/abs/2005ApJ...625..575S
    ForegroundComponent(130.0u"MHz", 0.014u"mK^2", 1.0, 2.10, 35)
end

function galactic_synchrotron()
    # http://adsabs.harvard.edu/abs/2005ApJ...625..575S
    ForegroundComponent(130.0u"MHz", 700u"mK^2", 2.4, 2.80, 4.0)
end

function galactic_free_free()
    # http://adsabs.harvard.edu/abs/2005ApJ...625..575S
    ForegroundComponent(130.0u"MHz", 0.088"mK^2", 3.0, 2.15, 35)
end

