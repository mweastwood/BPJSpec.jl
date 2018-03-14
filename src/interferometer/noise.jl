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
    struct NoiseModel

This type represents the thermal noise contributed to the measurement of a set of $m$-modes.

A careful reading of Taylor, Carilli, Perley chapter 9 reveals that under the convention that
Stokes-I is $(\rm xx + yy)/2$ we get the following expressions:

Uncertainty on single polarization visibilities (in flux density units):
```math
σ_{\rm xx} = \frac{\sqrt{2} k T_{\rm sys}}{A_e sqrt(Δν\,τ)}
```

Uncertainty on Stokes I visibilities (in flux density units):
```math
σ_{\rm Stokes-I} = \frac{σ_{\rm xx}}{\sqrt{2}} = \frac{k T_{\rm sys}}{A_e sqrt(Δν\,τ)}
```

Where $k$ is the Boltzmann constant, $T_{\rm sys}$ is the system temperature, $A_e$ is the effective
collecting area, $Δν$ is the bandwidth, and $τ$ is the integration time.

However, for a dipole antenna, the effective collecting area is not a very physically meaningful
value. However, it turns out that we can relate the effective collecting are to the solid angle
subtended by the primary beam $Ω$:
```math
A_e = \frac{λ^2}{Ω}
```

!!! note
    There seems to be some ambiguity in the literature in regards to notation. I believe we
    originally assumed that $A_e$ refers to the maximum effective collecting area, and that we have
    normalized the beam to be unity in that direction.

Finally we end  up with the following expression after including an additional contribution due to
time smearing:
```math
σ_{\rm m-modes} = \frac{k T_{\rm sys} Ω}{λ^2 \sqrt{Δν\,τ\,N_{\rm int}}}
                  \sinc\left(\frac{m\tau}{\text{sidereal day}}\right)
```

**Fields:**

* `Tsys` specifies the system temperature
* `τ` specfies the length of a single integration
* `Nint` specifies the total number of integrations used in the dataset
* `Ω` is the solid angle subtended by the primary beam

**Usage:**

```jldoctest
julia> model = BPJSpec.NoiseModel(1000u"K", 13u"s", 6628, 2.41u"sr")
NoiseModel(Tsys = 1000.0 K, τ = 13.0 s, Nint = 6628, Ω = 2.410 sr)

julia> model(100, 74u"MHz", 24u"kHz")
4.456470453155544 Jy
```
"""
struct NoiseModel
    # TODO: the system temperature should increase at lower frequencies
    Tsys  :: typeof(1.0u"K")  # system temperature
    τ     :: typeof(1.0u"s")  # integration time (for a single time slice)
    Nint  :: Int              # total number of integrations used in the dataset
    Ω     :: typeof(1.0u"sr") # the solid angle of the primary beam
end

function Base.show(io::IO, model::NoiseModel)
    @printf(io, "NoiseModel(Tsys = %.1f K, τ = %.1f s, Nint = %d, Ω = %.3f sr)",
            u(u"K", model.Tsys), u(u"s", model.τ), model.Nint, u(u"sr", model.Ω))
end

"Compute the standard error of a visibility from the system temperature."
function standard_error(Tsys, ν, Δν, τ, Nint, Ω)
    λ  = uconvert(u"m", u"c"/ν)       # wavelength
    Ae = uconvert(u"m^2", λ^2/Ω)      # effective collecting area
    N  = uconvert(NoUnits, Δν*τ*Nint) # number of independent samples
    σ  = u"k"*Tsys/(Ae*√N)            # standard error of a visibility
    uconvert(u"Jy", σ)
end

"Compute the effect of time smearing on the noise covariance of m-modes (Shaw 2015)."
function time_smearing(m, τ)
    t_sidereal = 86164.09054u"s" # one sidereal day
    x = uconvert(NoUnits, m*τ/t_sidereal)
    # N.B. Julia defines sinc(x) = sin(πx)/(πx)
    sinc(x)
end

function (model::NoiseModel)(m, ν, Δν)
    σ = standard_error(model.Tsys, ν, Δν, model.τ, model.Nint, model.Ω)
    uconvert(u"Jy", σ * time_smearing(m, model.τ))
end

