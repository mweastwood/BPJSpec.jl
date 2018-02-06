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
    Tsys  :: typeof(1.0*u"K")   # system temperature
    Δν    :: typeof(1.0*u"kHz") # bandwidth
    τ     :: typeof(1.0*u"s")   # integration time (for a single time slice)
    Nint  :: Int                # total number of integrations used in the dataset
    Nbase :: Int                # number of baselines
end

struct NoiseMatrix <: BlockMatrix
    path        :: String
    progressbar :: Bool
    distribute  :: Bool
    cached      :: Ref{Bool}
    mmax        :: Int
    frequencies :: Vector{typeof(1.0*u"Hz")}
    blocks      :: Matrix{Diagonal{Float64}}
    model       :: NoiseModel

    function NoiseMatrix(path, mmax, frequencies, model, write=true;
                         progressbar=false, distribute=false, cached=false)
        if write
            isdir(path) || mkpath(path)
            save(joinpath(path, "METADATA.jld2"), "mmax", mmax, "frequencies", frequencies)
        end
        blocks = Array{Matrix{Complex128}}(mmax+1, length(frequencies))
        output = new(path, progressbar, distribute, Ref(cached), mmax, frequencies, blocks, model)
        if write
            compute!(output)
        elseif cached
            cache!(output)
        end
        output
    end
end

function NoiseMatrix(path; kwargs...)
    mmax, frequencies = load(joinpath(path, "METADATA.jld2"), "mmax", "frequencies")
    NoiseMatrix(path, mmax, frequencies, false; kwargs...)
end

function compute!(matrix::NoiseMatrix)
    Nfreq = length(matrix.frequencies)
    model = matrix.model
    for β = 1:Nfreq
        ν = matrix.frequencies[β]
        σ = standard_error(model.Tsys, ν, model.Δν, model.τ, model.Nint)
        for m = 0:matrix.mmax
            σm = σ * time_smearing(m, model.τ)
            Nm = σm^2 .* ones(two(m)*model.Nbase)
            matrix[m, β] = Diagonal(Nm)
        end
    end
end

"Compute the standard error of a visibility from the system temperature."
function standard_error(Tsys, ν, Δν, τ, Nint)
    λ  = uconvert(u"m", u"c"/ν)       # wavelength
    Ae = uconvert(u"m^2", λ^2/(4π))   # effective collecting area (0th order approximation)
    N  = uconvert(NoUnits, Δν*τ*Nint) # number of independent samples
    σ  = u"k"*Tsys/(Ae*√N)            # standard error of a visibility
    ustrip(uconvert(u"Jy", σ))
end

"Compute the effect of time smearing on the noise covariance of m-modes (Shaw 2015)."
function time_smearing(m, τ)
    t_sidereal = 86164.09054u"s" # one sidereal day
    x = uconvert(NoUnits, m*τ/t_sidereal)
    # N.B. Julia defines sinc(x) = sin(πx)/(πx)
    sinc(x)
end

Base.show(io::IO, matrix::NoiseMatrix) = print(io, "NoiseMatrix: ", matrix.path)

indices(matrix::NoiseMatrix) =
    [(m, β) for β = 1:length(matrix.frequencies) for m = 0:matrix.mmax]

function Base.getindex(matrix::NoiseMatrix, m, β)
    if matrix.cached[]
        return matrix.blocks[m+1, β]
    else
        return read_from_disk(matrix, m, β)
    end
end

function Base.setindex!(matrix::NoiseMatrix, block, m, β)
    if matrix.cached[]
        matrix.blocks[m+1, β] = block
    else
        write_to_disk(matrix, block, m, β)
    end
    block
end

function cache!(matrix::NoiseMatrix)
    matrix.cached[] = true
    for β = 1:length(matrix.frequencies), m = 0:matrix.mmax
        matrix.blocks[m+1, β] = read_from_disk(matrix, m, β)
    end
    matrix
end

function flush!(matrix::NoiseMatrix)
    for β = 1:length(matrix.frequencies), m = 0:matrix.mmax
        write_to_disk(matrix, matrix.blocks[m+1, β], m, β)
    end
    matrix.cached[] = false
    matrix
end

function read_from_disk(matrix::NoiseMatrix, m, β)
    ν = matrix.frequencies[β]
    filename   = @sprintf("%.3fMHz.jld2", ustrip(uconvert(u"MHz", ν)))
    objectname = @sprintf("%04d", m)
    load(joinpath(matrix.path, filename), objectname) :: Diagonal{Float64}
end

function write_to_disk(matrix::NoiseMatrix, block::Diagonal{Float64}, m, β)
    ν = matrix.frequencies[β]
    filename   = @sprintf("%.3fMHz.jld2", ustrip(uconvert(u"MHz", ν)))
    objectname = @sprintf("%04d", m)
    jldopen(joinpath(matrix.path, filename), "a+") do file
        file[objectname] = block
    end
    block
end

