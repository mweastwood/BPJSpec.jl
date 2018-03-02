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

abstract type PowerSpectrum <: SkyComponent end

struct CylindricalPS <: PowerSpectrum
    zrange :: Tuple{Float64, Float64} # redshift range for this power spectrum
    kpara  :: Vector{typeof(1.0u"Mpc^-1")}
    kperp  :: Vector{typeof(1.0u"Mpc^-1")}
    power  :: Matrix{typeof(1.0u"K^2*Mpc^3")}
end

struct SphericalPS <: PowerSpectrum
    zrange :: Tuple{Float64, Float64} # redshift range for this power spectrum
    k      :: Vector{typeof(1.0u"Mpc^-1")}
    power  :: Vector{typeof(1.0u"K^2*Mpc^3")}
end

# It turns out that one of the slowest steps in calculating Cl(ν1, ν2) for a given power spectrum is
# simply computing the comoving radial distance to a given redshift. This function is actually
# pretty smooth so it can be well approximated by polynomials that can be rapidly evaluated.
function precomputation(model::PowerSpectrum)
    approximate(comoving_distance, model.zrange[1], model.zrange[2])
end

(model::CylindricalPS)(l, ν1, ν2) = Csignal(l, ν1, ν2, model)
(model::CylindricalPS)(l, ν1, ν2, comoving_distance_func) =
    Csignal(comoving_distance_func, l, ν1, ν2, model)

(model::SphericalPS)(l, ν1, ν2) = Csignal(l, ν1, ν2, model)
(model::SphericalPS)(l, ν1, ν2, comoving_distance_func) =
    Csignal(comoving_distance_func, l, ν1, ν2, model)


function Csignal(l, ν1, ν2, model)
    Csignal(comoving_distance, l, ν1, ν2, model)
end

function Csignal(comoving_distance_func, l, ν1, ν2, model)
    if ν1 == ν2
        return Csignal_same_redshift(comoving_distance_func, l, ν1, model)
    else
        return Csignal_different_redshift(comoving_distance_func, l, ν1, ν2, model)
    end
end

# The exact integral that needs to be performed is
#
#     Cₗ(ν₁, ν₂) = 2/π ∫ jₗ(k r₁) jₗ(k r₂) P(k) k² dk
#
# For the given frequencies ν₁ and ν₂ and corresponding comoving distances r₁ and r₂. The problem is
# that the spherical Bessel functions jₗ are highly oscillatory and this integral is very difficult
# to evaluate numerically.
#
# It turns out that the two spherical Bessel functions oscillate at slightly different frequencies
# and so there is a slow beat pattern imposed on the rapid oscillations. We can therefore
# approximate this integral as
#
#     Cₗ(ν₁, ν₂) = 1/(π r₁ r₂) ∫ cos(k Δr) P(k) k/sqrt(k² - k⟂²) dk
#
# Where k⟂ is l/r and r is some mean comoving distance to the given redshift (this is an ok
# approximation because we've already assumed that the comoving distances aren't so different when
# we right down a single power spectrum for both redshifts).
#
# But so far we've assumed a spherical power spectrum parameterized by the wave number k. In
# practice all the cool kids are doing cylindrically averaged power spectra parameterized by k⟂ and
# k∥. How do we generalize this integral in that case?
#
#     Cₗ(ν₁, ν₂) = 1/(π r₁ r₂) ∫ cos(k∥ Δr) P(k⟂, k∥) dk∥
#
# Ref: http://adsabs.harvard.edu/abs/2015PhRvD..91h3514S
#
# Now we're going to be given the power spectrum sampled on some grid. We will model the power
# spectrum as a piecewise linear interpolation of this grid (although we may want to look into a
# quadratic interpolation at some point).

function trapezoidal_rule_same_redshift(getx, gety, range, χ)
    out = 0.0u"K^2*Mpc^2"
    for i in range
        x1 = getx(i); x2 = getx(i+1)
        y1 = gety(i); y2 = gety(i+1)
        out += 0.5*(y1+y2)*(x2-x1)
    end
    uconvert(u"K^2", out/(π*χ^2))
end

function trapezoidal_rule_different_redshift(getx, gety, range, χ1, χ2)
    Δχ = χ2-χ1
    out = 0.0u"K^2*Mpc^2"
    for i in range
        x1 = getx(i); x2 = getx(i+1)
        y1 = gety(i); y2 = gety(i+1)
        out += (y2*sin(x2*Δχ) - y1*sin(x1*Δχ)) / Δχ
        out += (y2-y1) * (cos(x2*Δχ)-cos(x1*Δχ)) / (x2-x1) / Δχ^2
    end
    uconvert(u"K^2", out/(π*χ1*χ2))
end

function Csignal_same_redshift(comoving_distance_func, l, ν, model::CylindricalPS)
    z = redshift(ν)
    χ = comoving_distance_func(z)
    j, w1, w2 = interpolation_index_and_weights(model.kperp, l/χ)
    range   = 1:length(model.kpara)-1
    getx(i) = model.kpara[i]
    gety(i) = w1*model.power[i, j] + w2*model.power[i, j+1]
    trapezoidal_rule_same_redshift(getx, gety, range, χ)
end

function Csignal_different_redshift(comoving_distance_func, l, ν1, ν2, model::CylindricalPS)
    z1 = redshift(ν1)
    z2 = redshift(ν2)
    χ1 = comoving_distance_func(z1)
    χ2 = comoving_distance_func(z2)
    χ  = 0.5*(χ1+χ2)
    j, w1, w2 = interpolation_index_and_weights(model.kperp, l/χ)
    range   = 1:length(model.kpara)-1
    getx(i) = model.kpara[i]
    gety(i) = w1*model.power[i, j] + w2*model.power[i, j+1]
    trapezoidal_rule_different_redshift(getx, gety, range, χ1, χ2)
end

function Csignal_same_redshift(comoving_distance_func, l, ν, model::SphericalPS)
    z = redshift(ν)
    χ = comoving_distance_func(z)
    kperp = l/χ
    j, w1, w2 = interpolation_index_and_weights(model.k, kperp)
    start_power = w1*model.power[j] + w2*model.power[j+1]
    range   = j:length(model.k)-1
    getx(i) = i == j ? 0.0u"Mpc^-1" : sqrt(model.k[i]^2 - kperp^2)
    gety(i) = i == j ? start_power  : model.power[i]
    trapezoidal_rule_same_redshift(getx, gety, range, χ)
end

function Csignal_different_redshift(comoving_distance_func, l, ν1, ν2, model::SphericalPS)
    z1 = redshift(ν1)
    z2 = redshift(ν2)
    χ1 = comoving_distance_func(z1)
    χ2 = comoving_distance_func(z2)
    χ  = 0.5*(χ1+χ2)
    kperp = l/χ
    j, w1, w2 = interpolation_index_and_weights(model.k, kperp)
    start_power = w1*model.power[j] + w2*model.power[j+1]
    range   = j:length(model.k)-1
    getx(i) = i == j ? 0.0u"Mpc^-1" : sqrt(model.k[i]^2 - kperp^2)
    gety(i) = i == j ? start_power  : model.power[i]
    trapezoidal_rule_different_redshift(getx, gety, range, χ1, χ2)
end

function interpolation_index_and_weights(list, value)
    index = searchsortedlast(list, value)
    weight1 = (list[index+1]-value) / (list[index+1]-list[index])
    weight2 = 1 - weight1
    index, weight1, weight2
end

function fiducial_signal_model()
    kpara = logspace(log10(0.01), log10(1.0), 200) .* u"Mpc^-1"
    kperp = logspace(log10(0.01), log10(1.0), 200) .* u"Mpc^-1"
    unshift!(kpara, 0u"Mpc^-1")
    unshift!(kperp, 0u"Mpc^-1")
    k = sqrt.(kpara.^2 .+ kperp.'.^2)
    power = 2π^2 * 1000u"mK^2" ./ (k+0.1u"Mpc^-1").^3
    SignalModel((10., 30.), kpara, kperp, power)
end

