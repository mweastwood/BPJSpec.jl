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

abstract type Mechanism end

struct NoFile <: Mechanism end
Base.show(io::IO, ::NoFile) = print(io, "<no file>")
distribute_write(::NoFile) = false
distribute_read(::NoFile) = false

struct SingleFile <: Mechanism
    path :: String
end
Base.show(io::IO, storage::SingleFile) = print(io, storage.path)
magic(::Type{SingleFile}) = 0x0e85f27ea04cbee8 # rand(UInt64)
distribute_write(::SingleFile) = false
distribute_read(::SingleFile) = true

struct MultipleFiles <: Mechanism
    path :: String
end
Base.show(io::IO, storage::MultipleFiles) = print(io, storage.path)
magic(::Type{MultipleFiles}) = 0x24ef18b394f51d71 # rand(UInt64)
distribute_write(::MultipleFiles) = true
distribute_read(::MultipleFiles) = true

function write_metadata(storage::Mechanism, metadata)
    isdir(storage.path) || mkpath(storage.path)
    save(joinpath(storage.path, "METADATA.jld2"),
         "magic", magic(typeof(storage)), "metadata", metadata)
end

function read_metadata(path::String)
    magic′, metadata = load(joinpath(path, "METADATA.jld2"), "magic", "metadata")
    if     magic′ == magic(SingleFile)
        return SingleFile(path), metadata
    elseif magic′ == magic(MultipleFiles)
        return MultipleFiles(path), metadata
    else
        error("unknown magic number")
    end
end

function Base.setindex!(storage::SingleFile, block, idx::Int)
    filename   = "BLOCKS.jld2"
    objectname = @sprintf("%04d", idx)
    jldopen(joinpath(storage.path, filename), "a+") do file
        file[objectname] = block
    end
    block
end

function Base.getindex(storage::SingleFile, idx::Int)
    filename   = "BLOCKS.jld2"
    objectname = @sprintf("%04d", idx)
    load(joinpath(storage.path, filename), objectname)
end

function Base.setindex!(storage::SingleFile, block, idx::Int, jdx::Int)
    filename   = "BLOCKS.jld2"
    objectname = @sprintf("%04d/%04d", jdx, idx)
    jldopen(joinpath(storage.path, filename), "a+") do file
        file[objectname] = block
    end
    block
end

function Base.getindex(storage::SingleFile, idx::Int, jdx::Int)
    filename   = "BLOCKS.jld2"
    objectname = @sprintf("%04d/%04d", jdx, idx)
    load(joinpath(storage.path, filename), objectname)
end

function Base.setindex!(storage::MultipleFiles, block, idx::Int)
    filename   = @sprintf("%04d.jld2", idx)
    objectname = "block"
    save(joinpath(storage.path, filename), objectname, block)
    block
end

function Base.getindex(storage::MultipleFiles, idx::Int)
    filename   = @sprintf("%04d.jld2", idx)
    objectname = "block"
    load(joinpath(storage.path, filename), objectname)
end

function Base.setindex!(storage::MultipleFiles, block, idx::Int, jdx::Int)
    dirname    = @sprintf("%04d",      jdx)
    filename   = @sprintf("%04d.jld2", idx)
    objectname = "block"
    isdir(joinpath(storage.path, dirname)) || mkpath(joinpath(storage.path, dirname))
    save(joinpath(storage.path, dirname, filename), objectname, block)
    block
end

function Base.getindex(storage::MultipleFiles, idx::Int, jdx::Int)
    dirname    = @sprintf("%04d",      jdx)
    filename   = @sprintf("%04d.jld2", idx)
    objectname = "block"
    load(joinpath(storage.path, dirname, filename), objectname)
end

@inline Base.getindex(storage::Mechanism, idx::Tuple) = storage[idx...]
@inline Base.setindex!(storage::Mechanism, block, idx::Tuple) = storage[idx...] = block

struct Cache{T, N}
    used  :: Ref{Bool}
    cache :: Array{T, N}
end

function Cache{T}(length::Int) where T
    Cache{T, 1}(Ref(false), Vector{T}(length))
end

function Cache{T}(x::Int, y::Int) where T
    Cache{T, 2}(Ref(false), Matrix{T}(x, y))
end

function Cache{T}(tuple::Tuple) where T
    Cache{T}(tuple...)
end

Base.getindex(cache::Cache, idx) = cache.cache[idx]
Base.setindex!(cache::Cache, X, idx) = cache.cache[idx] = X
Base.getindex(cache::Cache, idx, jdx) = cache.cache[idx, jdx]
Base.setindex!(cache::Cache, X, idx, jdx) = cache.cache[idx, jdx] = X
@inline Base.getindex(cache::Cache, idx::Tuple) = cache[idx...]
@inline Base.setindex!(cache::Cache, X, idx::Tuple) = cache[idx...] = X
used(cache::Cache) = cache.used[]
set!(cache::Cache) = (cache.used[] = true; cache)
unset!(cache::Cache) = (cache.used[] = false; cache)

