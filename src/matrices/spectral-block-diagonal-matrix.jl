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
    Nfreq :: Int
    function SpectralBlockDiagonalMatrix(path, mmax, Nfreq)
        isdir(path) || mkpath(path)
        save(joinpath(path, "METADATA.jld2"), "mmax", mmax, "Nfreq", Nfreq)
        new(path, metadata)
    end
end

function SpectralBlockDiagonalMatrix(path)
    mmax, Nfreq = load(joinpath(path, "METADATA.jld2"), "mmax", "Nfreq")
    SpectralBlockDiagonalMatrix(path, mmax, Nfreq)
end

function Base.getindex(matrix::SpectralBlockDiagonalMatrix, m, channel)
    filename   = @sprintf("m=%04d.jld2", m)
    objectname = @sprintf("%04d", channel)
    load(joinpath(transfermatrix.path, filename), objectname) :: Matrix{Complex128}
end

function Base.setindex!(matrix::SpectralBlockDiagonalMatrix, block, m, channel)
    filename   = @sprintf("m=%04d.jld2", m)
    objectname = @sprintf("%04d", channel)
    jldopen(joinpath(transfermatrix.path, filename), "a+") do file
        file[objectname] = block
    end
    block
end

