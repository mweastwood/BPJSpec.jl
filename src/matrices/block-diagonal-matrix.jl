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

struct BlockDiagonalMatrix <: BlockMatrix
    path :: String
    mmax :: Int
    function BlockDiagonalMatrix(path, mmax, write=true)
        if write
            isdir(path) || mkpath(path)
            save(joinpath(path, "METADATA.jld2"), "mmax", mmax)
        end
        new(path, mmax)
    end
end

Base.show(io::IO, matrix::BlockDiagonalMatrix) = print(io, "BlockDiagonalMatrix: ", matrix.path)
Base.indices(matrix::BlockDiagonalMatrix) = 0:matrix.mmax

function BlockDiagonalMatrix(path)
    mmax = load(joinpath(path, "METADATA.jld2"), "mmax")
    BlockDiagonalMatrix(path, mmax, false)
end

function Base.getindex(matrix::BlockDiagonalMatrix, m)
    filename   = @sprintf("m=%4d.jld2", m)
    objectname = "block"
    load(joinpath(transfermatrix.path, filename), objectname) :: Matrix{Complex128}
end

function Base.setindex!(matrix::BlockDiagonalMatrix, block::Matrix{Complex128}, m)
    filename   = @sprintf("m=%4d.jld2", m)
    objectname = "block"
    save(joinpath(transfermatrix.path, filename), objectname, block)
    block
end

