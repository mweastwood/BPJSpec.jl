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

doc"""
    Cnoise(Tsys0,α,ν,ν0,Δν,τ_total,τ_int,m) -> Float64

Compute the expected variance of the $m$-modes due to thermal noise.

The system temperature is modeled as
\\[
    T_{sys} = T_{sys,0} \left(\frac{\nu}{\nu_0}\right)^{-\alpha},
\\]
$\Delta\nu$ gives the bandwidth, $\tau_{total}$ gives the total integration
time (in seconds), and $\tau_{int}$ gives the integration time corresponding
to a single integration (also in seconds).
"""
function Cnoise(Tsys0,α,ν,ν0,Δν,τ_total,τ_int,m)
    Tsys = Tsys0 * (ν/ν0)^(-α)
    tsid = 86164.09054 # sidereal day in seconds
    Tsys*Tsys / (τ_total*Δν) * sinc(m*τ_int/tsid)^2
end

immutable NoiseModel
    Tsys0::Float64
    α::Float64
    ν0::Float64
    Δν::Float64
    τ_total::Float64
    τ_int::Float64
end

function call(model::NoiseModel,m,ν)
    Cnoise(model.Tsys0,model.α,ν,model.ν0,model.Δν,model.τ_total,model.τ_int,m)
end

abstract SkyComponent

"""
    Cforeground(l,ν1,ν2,ν0,A,α,β,ζ) -> Float64

Evaluate a model for the multifrequency angular power spectrum of
a foreground component. Note that `ν0`, `A`, `α`, `β`, and `ζ` are
all parameters of the model.
"""
function Cforeground(l,ν1,ν2,ν0,A,α,β,ζ)
    (A * (l+1)^(-α)
        * (ν1*ν2/ν0^2)^(-β)
        * exp(-log(ν1/ν2)^2/(2*ζ^2)))
end

immutable ForegroundModel <: SkyComponent
    ν0::Float64
    A::Float64
    α::Float64
    β::Float64
    ζ::Float64
end

function call(model::ForegroundModel,l,ν1,ν2)
    Cforeground(l,ν1,ν2,model.ν0,model.A,model.α,model.β,model.ζ)
end

function Csignal(l,ν1,ν2,kpara,kperp,power)
    if ν1 == ν2
        z = redshift(ν1)
        χ = comoving_distance(z)

        # calculate weights for interpolating the value of the
        # power spectrum at k_perp = l/χ
        j = searchsortedlast(kperp,l/χ)
        weight_left  = 1 - (l/χ  -  kperp[j]) / (kperp[j+1] - kperp[j])
        weight_right = 1 - (kperp[j+1] - l/χ) / (kperp[j+1] - kperp[j])

        # integrate over k_para assuming that the power spectrum
        # is linear between the grid points
        out = 0.0
        for i = 1:length(kpara)-1
            k_start = kpara[i]
            k_stop  = kpara[i+1]
            P_start = weight_left*power[i,  j] + weight_right*power[i,  j+1]
            P_stop  = weight_left*power[i+1,j] + weight_right*power[i+1,j+1]
            out += 0.5*(P_start+P_stop)*(k_stop-k_start)
        end
        out /= π*χ^2
        return out
    else
        z1 = redshift(ν1)
        z2 = redshift(ν2)
        χ1 = comoving_distance(z1)
        χ2 = comoving_distance(z2)
        Δχ = χ2-χ1
        χ  = 0.5*(χ1+χ2)

        # calculate weights for interpolating the value of the
        # power spectrum at k_perp = l/χ
        j = searchsortedlast(kperp,l/χ)
        weight_left  = 1 - (l/χ  -  kperp[j]) / (kperp[j+1] - kperp[j])
        weight_right = 1 - (kperp[j+1] - l/χ) / (kperp[j+1] - kperp[j])

        # integrate over k_para assuming that the power spectrum
        # is linear between the grid points
        out = 0.0
        for i = 1:length(kpara)-1
            k_start = kpara[i]
            k_stop  = kpara[i+1]
            P_start = weight_left*power[i,  j] + weight_right*power[i,  j+1]
            P_stop  = weight_left*power[i+1,j] + weight_right*power[i+1,j+1]
            out += (P_stop*sin(k_stop*Δχ) - P_start*sin(k_start*Δχ)) / Δχ
            out += (P_stop-P_start) * (cos(k_stop*Δχ)-cos(k_start*Δχ)) / (k_stop-k_start) / Δχ^2
        end
        out /= π*χ1*χ2
        return out
    end
end

immutable SignalModel <: SkyComponent
    kpara::Vector{Float64} # Mpc⁻¹
    kperp::Vector{Float64} # Mpc⁻¹
    power::Matrix{Float64} # K² Mpc³
end

function dimensionful_powerspectrum(kpara,kperp,Δ²)
    out = similar(Δ²)
    for j = 1:length(kperp), i = 1:length(kpara)
        k = hypot(kpara[i],kperp[j])
        if k < eps(Float64)
            out[i,j] = 0
        else
            out[i,j] = 2π^2 * Δ²[i,j] / k^3
        end
    end
    out
end

function call(model::SignalModel,l,ν1,ν2)
    Csignal(l,ν1,ν2,model.kpara,model.kperp,model.power)
end

"""
    covariance_matrix(component::SkyComponent, ν, lmax, m)

Construct a covariance matrix for the given component of the sky.
"""
function covariance_matrix(component::SkyComponent, ν, lmax, m)
    Nfreq = length(ν)
    Ncoef = lmax-m+1 # number of spherical harmonic coefficients

    C = zeros(Complex128,Nfreq*Ncoef,Nfreq*Ncoef)
    for β2 = 1:Nfreq, β1 = 1:Nfreq
        block = zeros(Complex128,Ncoef,Ncoef)
        for l = m:lmax
            block[l-m+1,l-m+1] = component(l,ν[β1],ν[β2])
        end
        C[(β1-1)*Ncoef+1:β1*Ncoef,(β2-1)*Ncoef+1:β2*Ncoef] = block
    end
    C
end

