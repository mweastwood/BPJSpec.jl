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

function Csignal(l, ν1, ν2, kpara, kperp, power)
    Csignal(comoving_distance, l, ν1, ν2, kpara, kperp, power)
end

function Csignal(comoving_distance_func, l, ν1, ν2, kpara, kperp, power)
    if ν1 == ν2
        return Csignal_same_redshift(comoving_distance_func, l, ν1, kpara, kperp, power)
    else
        return Csignal_different_redshift(comoving_distance_func, l, ν1, ν2, kpara, kperp, power)
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

function Csignal_same_redshift(comoving_distance_func, l, ν, kpara, kperp, power)
    z = redshift(ν)
    χ = comoving_distance_func(z)

    # calculate weights for interpolating the value of the power spectrum at kperp = l/χ
    j = searchsortedlast(kperp, l/χ)
    weight_left  = 1 - (l/χ  -  kperp[j]) / (kperp[j+1] - kperp[j])
    weight_right = 1 - (kperp[j+1] - l/χ) / (kperp[j+1] - kperp[j])

    # integrate over kpara assuming that the power spectrum is linear between the grid points
    out = 0.0u"K^2*Mpc^2"
    for i = 1:length(kpara)-1
        k_start = kpara[i]
        k_stop  = kpara[i+1]
        P_start = weight_left*power[i,   j] + weight_right*power[i,   j+1]
        P_stop  = weight_left*power[i+1, j] + weight_right*power[i+1, j+1]
        out += 0.5*(P_start+P_stop)*(k_stop-k_start)
    end
    ustrip(uconvert(u"K^2", out/(π*χ^2)))
end

function Csignal_different_redshift(comoving_distance_func, l, ν1, ν2, kpara, kperp, power)
    z1 = redshift(ν1)
    z2 = redshift(ν2)
    χ1 = comoving_distance_func(z1)
    χ2 = comoving_distance_func(z2)
    Δχ = χ2-χ1
    χ  = 0.5*(χ1+χ2)

    # calculate weights for interpolating the value of the power spectrum at kperp = l/χ
    j = searchsortedlast(kperp, l/χ)
    weight_left  = 1 - (l/χ  -  kperp[j]) / (kperp[j+1] - kperp[j])
    weight_right = 1 - (kperp[j+1] - l/χ) / (kperp[j+1] - kperp[j])

    # integrate over kpara assuming that the power spectrum is linear between the grid points
    out = 0.0u"K^2*Mpc^2"
    for i = 1:length(kpara)-1
        k_start = kpara[i]
        k_stop  = kpara[i+1]
        P_start = weight_left*power[i,  j] + weight_right*power[i,  j+1]
        P_stop  = weight_left*power[i+1,j] + weight_right*power[i+1,j+1]
        out += (P_stop*sin(k_stop*Δχ) - P_start*sin(k_start*Δχ)) / Δχ
        out += (P_stop-P_start) * (cos(k_stop*Δχ)-cos(k_start*Δχ)) / (k_stop-k_start) / Δχ^2
    end
    ustrip(uconvert(u"K^2", out/(π*χ1*χ2)))
end

struct SignalModel <: SkyComponent
    zrange :: Tuple{Float64, Float64} # redshift range for this power spectrum
    kpara  :: Vector{typeof(1.0u"Mpc^-1")}
    kperp  :: Vector{typeof(1.0u"Mpc^-1")}
    power  :: Matrix{typeof(1.0u"K^2*Mpc^3")}
end

# It turns out that one of the slowest steps in calculating Cl(ν1, ν2) for a given power spectrum is
# simply computing the comoving radial distance to a given redshift. This function is actually
# pretty smooth so it can be well approximated by polynomials that can be rapidly evaluated.
function precomputation(model::SignalModel)
    approximate(comoving_distance, model.zrange[1], model.zrange[2])
end

function (model::SignalModel)(l, ν1, ν2)
    Csignal(l, ν1, ν2, model.kpara, model.kperp, model.power)
end

function (model::SignalModel)(l, ν1, ν2, comoving_distance_func)
    Csignal(comoving_distance_func, l, ν1, ν2, model.kpara, model.kperp, model.power)
end

function fiducial_signal_model()
    kpara = logspace(log10(0.01), log10(1.0), 200).*u"Mpc^-1"
    kperp = logspace(log10(0.01), log10(1.0), 200).*u"Mpc^-1"
    unshift!(kpara, 0u"Mpc^-1")
    unshift!(kperp, 0u"Mpc^-1")
    k = sqrt.(kpara.^2 .+ kperp.'.^2)
    power = 2π^2 * 1000u"mK^2" ./ (k+0.1u"Mpc^-1").^3
    SignalModel((10., 30.), kpara, kperp, power)
end

#function dimensionful_powerspectrum(kpara,kperp,Δ²)
#    out = similar(Δ²)
#    for j = 1:length(kperp), i = 1:length(kpara)
#        k = hypot(kpara[i],kperp[j])
#        if k < eps(Float64)
#            out[i,j] = 0
#        else
#            out[i,j] = 2π^2 * Δ²[i,j] / k^3
#        end
#    end
#    out
#end

