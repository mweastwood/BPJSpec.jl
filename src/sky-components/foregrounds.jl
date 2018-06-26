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
    struct ForegroundComponent <: SkyComponent

This type represents the contribution of a single foreground component to the multi-frequency
angular power spectrum.

```math
C_l(ν_1, ν_2) = A \left(\frac{l+1}{1001}\right)^{-α}
                  \left(\frac{ν_1 ν_2}{ν_0^2}\right)^{-β}
                  \exp\left(-\frac{\log^2(ν_1/ν_2)}{2ζ^2}\right)
```

**Fields:**

* `ν0` specifies the reference frequency
* `A` specifies the overall amplitude at the reference frequency
* `α` is the power-law index for the multipole moment $l$
* `β` is the power-law index for frequency
* `ζ` essentially determines how quickly the foreground component decorrelates in frequency

**Usage:**

Some quick foreground models can be constructed based on the work of [Santos, Cooray & Knox
2005](http://adsabs.harvard.edu/abs/2005ApJ...625..575S). These can be accessed with the functions
[`extragalactic_point_sources`](@ref), [`extragalactic_free_free`](@ref),
[`galactic_synchrotron`](@ref), and [`galactic_free_free`](@ref).

```jldoctest
julia> BPJSpec.extragalactic_point_sources()
ForegroundComponent(ν0 = 130.000 MHz, A = 57.000 mK², α = 1.100, β = 2.070, ζ = 1.000)

julia> BPJSpec.extragalactic_free_free()
ForegroundComponent(ν0 = 130.000 MHz, A = 0.014 mK², α = 1.000, β = 2.100, ζ = 35.000)

julia> BPJSpec.galactic_synchrotron()
ForegroundComponent(ν0 = 130.000 MHz, A = 700.000 mK², α = 2.400, β = 2.800, ζ = 4.000)

julia> BPJSpec.galactic_free_free()
ForegroundComponent(ν0 = 130.000 MHz, A = 0.088 mK², α = 3.000, β = 2.150, ζ = 35.000)

julia> component = BPJSpec.ForegroundComponent(100u"MHz", 1u"K^2", 1, 1, 100)
       component(100, 74u"MHz", 76u"MHz")
17.622494197510505 K^2
```
"""
struct ForegroundComponent <: SkyComponent
    ν0 :: typeof(1.0u"Hz")
    A  :: typeof(1.0u"K^2")
    α  :: Float64
    β  :: Float64
    ζ  :: Float64
end

function Base.show(io::IO, component::ForegroundComponent)
    @printf(io, "ForegroundComponent(ν0 = %.3f MHz, A = %.3f mK², α = %.3f, β = %.3f, ζ = %.3f)",
            u(u"MHz", component.ν0), u(u"mK^2", component.A), component.α, component.β, component.ζ)
end

function (component::ForegroundComponent)(l, ν1, ν2)
    Cforeground(l, ν1, ν2, component.ν0, component.A, component.α, component.β, component.ζ)
end

function Cforeground(l, ν1, ν2, ν0, A, α, β, ζ)
    x = (l+1)/1001
    y = uconvert(Unitful.NoUnits, ν1*ν2/ν0^2)
    z = log(ν1/ν2)
    uconvert(u"K^2", A * x^-α * y^-β * exp(-z^2/(2ζ^2)))
end

"""
    extragalactic_point_sources()

Constructs a model of extragalactic point sources based on the work of [Santos, Cooray & Knox
2005](http://adsabs.harvard.edu/abs/2005ApJ...625..575S).
"""
function extragalactic_point_sources()
    # http://adsabs.harvard.edu/abs/2005ApJ...625..575S
    ForegroundComponent(130.0u"MHz", 57.0u"mK^2", 1.1, 2.07, 1.0)
end

"""
    extragalactic_free_free()

Constructs a model of extragalactic free-free emission based on the work of [Santos, Cooray & Knox
2005](http://adsabs.harvard.edu/abs/2005ApJ...625..575S).
"""
function extragalactic_free_free()
    # http://adsabs.harvard.edu/abs/2005ApJ...625..575S
    ForegroundComponent(130.0u"MHz", 0.014u"mK^2", 1.0, 2.10, 35)
end

"""
    galactic_synchrotron()

Constructs a model of galactic synchrotron emission based on the work of [Santos, Cooray & Knox
2005](http://adsabs.harvard.edu/abs/2005ApJ...625..575S).
"""
function galactic_synchrotron()
    # http://adsabs.harvard.edu/abs/2005ApJ...625..575S
    ForegroundComponent(130.0u"MHz", 700u"mK^2", 2.4, 2.80, 4.0)
end

"""
    galactic_free_free()

Constructs a model of galactic free-free emission based on the work of [Santos, Cooray & Knox
2005](http://adsabs.harvard.edu/abs/2005ApJ...625..575S).
"""
function galactic_free_free()
    # http://adsabs.harvard.edu/abs/2005ApJ...625..575S
    ForegroundComponent(130.0u"MHz", 0.088u"mK^2", 3.0, 2.15, 35)
end

doc"""
    struct GeneralizedForegroundComponent <: SkyComponent

This type allows for a more general representation of the foreground emission. Instead of using a
parametric model, the multi-frequency angular power spectrum is discretized.

Furthermore, we assume that the power spectrum is a function of only the geometric mean of the two
frequencies (ie. $\sqrt{\nu_1 \nu_2}$). This is driven by the expectation that there is very little
decorrelation of the foreground emission between frequency channels. Instead, this allows for a
power-law spectrum of the foreground emission.
"""
struct GeneralForegroundComponent <: SkyComponent
    l :: Vector{Int}
    ν :: Vector{typeof(1.0u"Hz")}
    amplitude :: Matrix{typeof(1.0u"K^2")}
end

function (component::GeneralForegroundComponent)(l, ν1, ν2)
    ν = sqrt(ν1 * ν2) # geometric mean
    idx = searchsortedlast(component.l, l)
    jdx = searchsortedlast(component.ν, ν)
    weight_l = (component.l[idx+1] - l) / (component.l[idx+1] - component.l[idx])
    weight_ν = (component.ν[idx+1] - ν) / (component.ν[idx+1] - component.ν[idx])
    (component.amplitude[idx, jdx] * weight_l * weight_ν
        + component.amplitude[idx+1, jdx] * (1-weight_l) * weight_ν
        + component.amplitude[idx, jdx+1] * weight_l * (1-weight_ν)
        + component.amplitude[idx+1, jdx+1] * (1-weight_l) * (1-weight-ν))
end

