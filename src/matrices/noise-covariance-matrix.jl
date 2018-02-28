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

struct NoiseModel
    # TODO: the system temperature should increase at lower frequencies
    Tsys  :: typeof(1.0*u"K")  # system temperature
    τ     :: typeof(1.0*u"s")  # integration time (for a single time slice)
    Nint  :: Int               # total number of integrations used in the dataset
    Ω     :: typeof(1.0*u"sr") # the solid angle of the primary beam
end

const NoiseCovarianceMatrix = SpectralBlockDiagonalMatrix{Diagonal{Float64}}

function NoiseCovarianceMatrix(path, mmax, frequencies, bandwidth, hierarchy, noise::NoiseModel)
    output = NoiseCovarianceMatrix(path, mmax, frequencies, bandwidth)
    compute!(output, hierarchy, noise)
    output
end

function compute!(matrix::NoiseCovarianceMatrix, hierarchy, noise::NoiseModel)
    Nfreq = length(matrix.frequencies)
    for β = 1:Nfreq
        ν  = matrix.frequencies[β]
        Δν = matrix.bandwidth[β]
        σ  = ustrip(uconvert(u"Jy", standard_error(noise.Tsys, ν, Δν, noise.τ,
                                                   noise.Nint, noise.Ω)))
        for m = 0:matrix.mmax
            σm = σ * time_smearing(m, noise.τ)
            Nm = σm^2 .* ones(two(m)*Nbase(hierarchy, m))
            matrix[m, β] = Diagonal(Nm)
        end
    end
end

"Compute the standard error of a visibility from the system temperature."
function standard_error(Tsys, ν, Δν, τ, Nint, Ω)
    # A careful reading of Taylor, Carilli, Perley chapter 9 reveals that under the convention that
    # Stokes I is 0.5*(xx + yy) we get the following expressions:
    #
    # Uncertainty on single polarization visibilities (in flux density units):
    #     σₓₓ = sqrt(2) k T / A sqrt(Δν τ)
    #
    # Uncertainty on Stokes I visibilities (in flux density units):
    #     σ = σₓₓ / sqrt(2) = k T / A sqrt(Δν τ)
    #
    # Where k = Boltzman constant
    #       T = system temperature
    #       A = effective collecting area
    #       Δν = bandwidth
    #       τ = total integration time
    #
    # Now there's still a question of what do we mean by "effective collecting area"? There are two
    # possibilities here:
    # 1. The angle averaged effective collecting area = λ^2/(4π)
    # 2. The peak effective collecting area = λ^2/Ω where Ω is the beam solid angle
    #
    # I am of the opinion that we actually want option 2 here, but this is something I should
    # revisit and get a couple of other opinions (TODO!)
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

Base.show(io::IO, matrix::NoiseCovarianceMatrix) =
    print(io, "NoiseCovarianceMatrix: ", matrix.path)

