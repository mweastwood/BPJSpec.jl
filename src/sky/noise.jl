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

#immutable NoiseMeta <: Metadata
#    m::UnitRange{Int}
#    ν::Vector{Float64}
#
#    function NoiseMeta(m,ν)
#        if length(m) > 1 && length(ν) > 1
#            error("Cannot simultaneously have multiple values of m and multiple frequency channels.")
#        end
#        new(m,ν)
#    end
#end
#
#NoiseMeta(m::Int,ν::AbstractVector) = NoiseMeta(m:m,collect(ν))
#NoiseMeta(mmax::Int,ν::Float64) = NoiseMeta(0:mmax,[ν])
#
#==(lhs::NoiseMeta,rhs::NoiseMeta) = lhs.m == rhs.m && lhs.ν == rhs.ν
#
#typealias DiagonalNoiseMatrix Blocks{DiagonalMatrixBlock, NoiseMeta}
#typealias NoiseMatrix Blocks{MatrixBlock, NoiseMeta}
#
#initial_block_size(::Type{DiagonalNoiseMatrix}, Nbase, m) = (two(m)*Nbase,)
#
#function call(::Type{DiagonalNoiseMatrix}, Nbase::Int, mmax::Int, ν::Float64)
#    meta = NoiseMeta(mmax,ν)
#    blocks = DiagonalMatrixBlock[]
#    for m = 0:mmax
#        sz = initial_block_size(DiagonalNoiseMatrix,Nbase,m)
#        push!(blocks,DiagonalMatrixBlock(sz))
#    end
#    DiagonalNoiseMatrix(blocks,meta)
#end
#
#is_single_frequency(meta::NoiseMeta) = length(meta.ν) == 1
#is_single_m(meta::NoiseMeta) = length(meta.m) == 1
#is_single_frequency(N::NoiseMatrix) = is_single_frequency(N.meta)
#is_single_m(N::NoiseMatrix) = is_single_m(N.meta)
#
#mmax(meta::NoiseMeta) = maximum(meta.m)
#mmax(N::NoiseMatrix) = mmax(N.meta)
#
#Nfreq(meta::NoiseMeta) = length(meta.ν)
#Nfreq(N::NoiseMatrix) = length(N.meta)
#
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
#
#function call(model::NoiseModel,m,ν)
#    Cnoise(model.Tsys0,model.α,ν,model.ν0,model.Δν,model.τ_total,model.τ_int,m)
#end
#
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
#
#function save(filename, N::NoiseMatrix)
#    if !isfile(filename)
#        jldopen(filename,"w",compress=true) do file
#            file["description"] = "noise covariance matrix"
#        end
#    end
#
#    jldopen(filename,"r+",compress=true) do file
#        if read(file["description"]) != "noise covariance matrix"
#            error("Attempting to write to a file that does not contain m-modes.")
#        end
#
#        if is_single_frequency(N)
#            ν = N.meta.ν[1]
#            name = @sprintf("%.3fMHz",ν/1e6)
#            name in names(file) || g_create(file,name)
#            group = file[name]
#            for m = 0:mmax(N)
#                block = N[m+1]
#                group[string(m)] = block.block
#            end
#        elseif is_single_m(N)
#            m = N.meta.m[1]
#            for β = 1:Nfreq(N)
#                ν = N.meta.ν[β]
#                name = @sprintf("%.3fMHz",ν/1e6)
#                name in names(file) || g_create(file,name)
#                group = file[name]
#                block = N[β]
#                group[string(m)] = block.block
#            end
#        end
#    end
#end
#
#function load(filename, meta::NoiseMeta)
#    blocks = MatrixBlock[]
#    jldopen(filename,"r") do file
#        if read(file["description"]) != "noise covariance matrix"
#            error("Attempting to read from a file that does not contain m-modes.")
#        end
#
#        if is_single_frequency(meta)
#            ν = meta.ν[1]
#            name = @sprintf("%.3fMHz",ν/1e6)
#            group = file[name]
#            for m = 0:mmax(meta)
#                block = group[string(m)] |> read
#                push!(blocks,MatrixBlock(block))
#            end
#        elseif is_single_m(meta)
#            m = meta.m[1]
#            for β = 1:Nfreq(meta)
#                ν = meta.ν[β]
#                name = @sprintf("%.3fMHz",ν/1e6)
#                group = file[name]
#                block = group[string(m)] |> read
#                push!(blocks,MatrixBlock(block))
#            end
#        end
#    end
#    NoiseMatrix(blocks,meta)
#end

