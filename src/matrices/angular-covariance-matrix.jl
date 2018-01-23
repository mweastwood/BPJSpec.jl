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
    path :: String
    lmax :: Int
    function AngularCovarianceMatrix(path, lmax)
        isdir(path) || mkpath(path)
        save(joinpath(path, "METADATA.jld2"), "lmax", lmax)
        new(path, lmax)
    end
end

function AngularCovarianceMatrix(path)
    lmax = load(joinpath(path, "METADATA.jld2"), "lmax")
    AngularCovarianceMatrix(path, lmax)
end

function Base.getindex(matrix::AngularCovarianceMatrix, l)
    filename   = @sprintf("l=%4d.jld2", l)
    objectname = "block"
    load(joinpath(transfermatrix.path, filename), objectname) :: Matrix{Float64}
end

function Base.setindex!(matrix::AngularCovarianceMatrix, block, l)
    filename   = @sprintf("l=%4d.jld2", l)
    objectname = "block"
    save(joinpath(transfermatrix.path, filename), objectname, block)
    block
end

