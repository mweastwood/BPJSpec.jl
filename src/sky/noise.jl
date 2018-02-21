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

#doc"""
#    Cnoise(Tsys0,α,ν,ν0,Δν,τ_total,τ_int,m) -> Float64
#
#Compute the expected variance of the $m$-modes due to thermal noise.
#
#The system temperature is modeled as
#\\[
#    T_{sys} = T_{sys,0} \left(\frac{\nu}{\nu_0}\right)^{-\alpha},
#\\]
#$\Delta\nu$ gives the bandwidth, $\tau_{total}$ gives the total integration
#time (in seconds), and $\tau_{int}$ gives the integration time corresponding
#to a single integration (also in seconds).
#"""
#function Cnoise(Tsys0,α,ν,ν0,Δν,τ_total,τ_int,m)
#    Tsys = Tsys0 * (ν/ν0)^(-α)
#    tsid = 86164.09054 # sidereal day in seconds
#    Tsys*Tsys / (τ_total*Δν) * sinc(m*τ_int/tsid)^2
#end
#
#immutable NoiseModel
#    Tsys0::Float64
#    α::Float64
#    ν0::Float64
#    Δν::Float64
#    τ_total::Float64
#    τ_int::Float64
#end

#function covariance_matrix(noise::NoiseModel, v::MModes)
#    is_single_frequency(v) || error("Expected single-frequency m-modes.")
#    ν = v.meta.ν[1]
#    Nbase = length(v[1])
#    N = DiagonalNoiseMatrix(Nbase,mmax(v),ν)
#    for m = 0:mmax(v)
#        amplitude = noise(m,ν)
#        for α = 1:size(N[m+1])[1]
#            N[m+1][α] = amplitude
#        end
#    end
#    N
#end

