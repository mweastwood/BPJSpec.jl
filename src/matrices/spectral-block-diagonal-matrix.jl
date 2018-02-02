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

struct SpectralBlockDiagonalMatrix <: BlockMatrix
    path  :: String
    mmax  :: Int
    frequencies :: Vector{typeof(1.0*u"Hz")}
    function SpectralBlockDiagonalMatrix(path, mmax, frequencies, write=true)
        if write
            isdir(path) || mkpath(path)
            save(joinpath(path, "METADATA.jld2"), "mmax", mmax, "frequencies", frequencies)
        end
        new(path, mmax, frequencies)
    end
end

function SpectralBlockDiagonalMatrix(path)
    mmax, frequencies = load(joinpath(path, "METADATA.jld2"), "mmax", "frequencies")
    SpectralBlockDiagonalMatrix(path, mmax, frequencies, false)
end

Base.show(io::IO, matrix::SpectralBlockDiagonalMatrix) =
    print(io, "SpectralBlockDiagonalMatrix: ", matrix.path)
Base.indices(matrix::SpectralBlockDiagonalMatrix) = (0:matrix.mmax, matrix.frequencies)

function Base.getindex(matrix::SpectralBlockDiagonalMatrix, m, ν)
    if !(uconvert(u"Hz", ν) in matrix.frequencies)
        error("unkown frequency")
    end
    filename   = @sprintf("%.3fMHz.jld2", ustrip(uconvert(u"MHz", ν)))
    objectname = @sprintf("%04d", m)
    load(joinpath(matrix.path, filename), objectname) :: Matrix{Complex128}
end

function Base.setindex!(matrix::SpectralBlockDiagonalMatrix, block, m, ν)
    if !(uconvert(u"Hz", ν) in matrix.frequencies)
        error("unkown frequency")
    end
    filename   = @sprintf("%.3fMHz.jld2", ustrip(uconvert(u"MHz", ν)))
    objectname = @sprintf("%04d", m)
    jldopen(joinpath(matrix.path, filename), "a+") do file
        file[objectname] = block
    end
    block
end

function Base.getindex(matrix::SpectralBlockDiagonalMatrix, m)
    blocks = [matrix[m, ν] for ν in matrix.frequencies]
    X = sum(size.(blocks, 1))
    Y = sum(size.(blocks, 2))
    output = zeros(Complex128, X, Y)
    x = y = 1
    for block in blocks
        output[x:x+size(block, 1)-1, y:y+size(block, 2)-1] = block
        x += size(block, 1)
        y += size(block, 2)
    end
    output
end

