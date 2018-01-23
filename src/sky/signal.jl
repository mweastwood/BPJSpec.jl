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
    if ν1 == ν2
        return Csignal_same_redshift(l, ν1, kpara, kperp, power)
    else
        return Csignal_different_redshift(l, ν1, ν2, kpara, kperp, power)
    end
end

function Csignal_same_redshift(l, ν, kpara, kperp, power)
    z = redshift(ν)
    χ = comoving_distance(z)

    # calculate weights for interpolating the value of the power spectrum at k_perp = l/χ
    j = searchsortedlast(kperp,l/χ)
    weight_left  = 1 - (l/χ  -  kperp[j]) / (kperp[j+1] - kperp[j])
    weight_right = 1 - (kperp[j+1] - l/χ) / (kperp[j+1] - kperp[j])

    # integrate over k_para assuming that the power spectrum is linear between the grid points
    out = 0.0
    for i = 1:length(kpara)-1
        k_start = kpara[i]
        k_stop  = kpara[i+1]
        P_start = weight_left*power[i,  j] + weight_right*power[i,  j+1]
        P_stop  = weight_left*power[i+1,j] + weight_right*power[i+1,j+1]
        out += 0.5*(P_start+P_stop)*(k_stop-k_start)
    end
    out /= π*χ^2
    out
end

function Csignal_different_redshift(l, ν1, ν2, kpara, kperp, power)
    z1 = redshift(ν1)
    z2 = redshift(ν2)
    χ1 = comoving_distance(z1)
    χ2 = comoving_distance(z2)
    Δχ = χ2-χ1
    χ  = 0.5*(χ1+χ2)

    # calculate weights for interpolating the value of the power spectrum at k_perp = l/χ
    j = searchsortedlast(kperp,l/χ)
    weight_left  = 1 - (l/χ  -  kperp[j]) / (kperp[j+1] - kperp[j])
    weight_right = 1 - (kperp[j+1] - l/χ) / (kperp[j+1] - kperp[j])

    # integrate over k_para assuming that the power spectrum is linear between the grid points
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
    out
end

#immutable SignalModel <: SkyComponent
#    kpara :: Vector{Float64} # Mpc⁻¹
#    kperp :: Vector{Float64} # Mpc⁻¹
#    power :: Matrix{Float64} # K² Mpc³
#end
#
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
#
#function call(model::SignalModel,l,ν1,ν2)
#    Csignal(l,ν1,ν2,model.kpara,model.kperp,model.power)
#end

