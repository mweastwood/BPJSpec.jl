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
        new(path, lmax, frequencies, component)
    end
end

function AngularCovarianceMatrix(path)
    lmax, frequencies, component = load(joinpath(path, "METADATA.jld2"),
                                        "lmax", "frequencies", "component")
    AngularCovarianceMatrix(path, lmax, frequencies, component)
end

function Base.getindex(matrix::AngularCovarianceMatrix, l)
    filename   = @sprintf("l=%04d.jld2", l)
    objectname = "block"
    load(joinpath(matrix.path, filename), objectname) :: Matrix{Float64}
end

function Base.setindex!(matrix::AngularCovarianceMatrix, block, l)
    filename   = @sprintf("l=%04d.jld2", l)
    objectname = "block"
    save(joinpath(matrix.path, filename), objectname, block)
    block
end

function compute!(matrix::AngularCovarianceMatrix)
    Nfreq = length(matrix.frequencies)
    prg = Progress(matrix.lmax+1)
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
        next!(prg)
    end
end

function densify(matrix::AngularCovarianceMatrix, m)
    Nfreq = length(matrix.frequencies)
    lmax  = matrix.lmax
    N = (lmax-m+1)*Nfreq
    output = zeros(Float64, N, N)
    for l = m:lmax
        block = matrix[l]
        for β1 = 1:Nfreq, β2 = 1:Nfreq
            x = (lmax-m+1)*(β1-1) + l + 1
            y = (lmax-m+1)*(β2-1) + l + 1
            output[x, y] = block[β1, β2]
        end
    end
    output
end

