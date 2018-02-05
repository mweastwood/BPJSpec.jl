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

struct AngularCovarianceMatrix <: BlockMatrix
    path        :: String
    progressbar :: Bool
    distribute  :: Bool
    cached      :: Ref{Bool}
    lmax        :: Int
    frequencies :: Vector{typeof(1.0*u"Hz")}
    component   :: SkyComponent
    blocks      :: Vector{Matrix{Float64}}

    function AngularCovarianceMatrix(path, lmax, frequencies, component, write=true;
                                     progressbar=false, distribute=false, cached=false)
        if write
            isdir(path) || mkpath(path)
            save(joinpath(path, "METADATA.jld2"), "lmax", lmax,
                 "frequencies", frequencies, "component", component)
        end
        blocks = Matrix{Float64}[]
        output = new(path, progressbar, distribute, Ref(cached),
                     lmax, frequencies, component, blocks)
        if write
            compute!(output)
        elseif cached
            cache!(output)
        end
        output
    end
end

function AngularCovarianceMatrix(path; kwargs...)
    lmax, frequencies, component = load(joinpath(path, "METADATA.jld2"),
                                        "lmax", "frequencies", "component")
    AngularCovarianceMatrix(path, lmax, frequencies, component, false; kwargs...)
end

Base.show(io::IO, matrix::AngularCovarianceMatrix) =
    print(io, "AngularCovarianceMatrix: ", matrix.path)

indices(matrix::AngularCovarianceMatrix) =
    [(l, m) for m = 0:matrix.mmax for l = m:matrix.lmax]

function Base.getindex(matrix::AngularCovarianceMatrix, l::Integer)
    if matrix.cached[]
        return matrix.blocks[l+1]
    else
        return read_from_disk(matrix, l)
    end
end

function Base.setindex!(matrix::AngularCovarianceMatrix, block, l::Integer)
    if matrix.cached[]
        matrix.blocks[l+1] = block
    else
        write_to_disk(matrix, block, l)
    end
    block
end

Base.getindex(matrix::AngularCovarianceMatrix, l, m) = matrix[l]
Base.setindex!(matrix::AngularCovarianceMatrix, block, l, m) = matrix[l] = block

function compute!(matrix::AngularCovarianceMatrix)
    Nfreq = length(matrix.frequencies)
    for l = 0:matrix.lmax
        block = zeros(Float64, Nfreq, Nfreq)
        for β1 = 1:Nfreq
            ν1 = matrix.frequencies[β1]
            block[β1, β1] = matrix.component(l, ν1, ν1)
            for β2 = β1+1:Nfreq
                ν2 = matrix.frequencies[β2]
                block[β1, β2] = matrix.component(l, ν1, ν2)
                block[β2, β1] = block[β1, β2]
            end
        end
        matrix[l] = block
    end
end

function cache!(matrix::AngularCovarianceMatrix)
    matrix.cached[] = true
    empty!(matrix.blocks)
    for l = 0:matrix.lmax
        push!(matrix.blocks, read_from_disk(matrix, l))
    end
    matrix
end

function flush!(matrix::AngularCovarianceMatrix)
    for l = 0:matrix.lmax
        write_to_disk(matrix, matrix.blocks[l+1], l)
    end
    empty!(matrix.blocks)
    matrix.cached[] = false
    matrix
end

function read_from_disk(matrix::AngularCovarianceMatrix, l::Integer)
    filename   = @sprintf("l=%04d.jld2", l)
    objectname = "block"
    load(joinpath(matrix.path, filename), objectname) :: Matrix{Float64}
end

function write_to_disk(matrix::AngularCovarianceMatrix, block::Matrix{Float64}, l::Integer)
    filename   = @sprintf("l=%04d.jld2", l)
    objectname = "block"
    save(joinpath(matrix.path, filename), objectname, block)
    block
end

