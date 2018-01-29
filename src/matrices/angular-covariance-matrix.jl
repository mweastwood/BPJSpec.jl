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
    frequencies :: Vector{typeof(1.0*u"Hz")}
    component   :: SkyComponent
    function AngularCovarianceMatrix(path, lmax, frequencies, component)
        isdir(path) || mkpath(path)
        save(joinpath(path, "METADATA.jld2"), "lmax", lmax,
             "frequencies", frequencies, "component", component)
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
                block[β2, β1] = block[β2, β1]
            end
        end
        matrix[l] = block
    end
end

