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
NoFile(path) = NoFile()
Base.show(io::IO, ::NoFile) = print(io, "<no file>")
distribute_write(::NoFile) = false
distribute_read(::NoFile) = false

struct SingleFile <: Mechanism
    path :: String
end
Base.show(io::IO, storage::SingleFile) = print(io, "\"", storage.path, "\"")
distribute_write(::SingleFile) = false
distribute_read(::SingleFile) = true

struct MultipleFiles <: Mechanism
    path :: String
end
Base.show(io::IO, storage::MultipleFiles) = print(io, "\"", storage.path, "/*\"")
distribute_write(::MultipleFiles) = true
distribute_read(::MultipleFiles) = true

struct HierarchicalStorage <: Mechanism
    path :: String
    hierarchy :: Hierarchy
end
Base.show(io::IO, storage::HierarchicalStorage) = print(io, "\"", storage.path, "/+\"")
distribute_write(::HierarchicalStorage) = false
distribute_read(::HierarchicalStorage) = true

function write_metadata(storage::Mechanism, metadata)
    isdir(storage.path) || mkpath(storage.path)
    jldopen(joinpath(storage.path, "METADATA.jld2"), mode[w]..., IOStream) do file
        file["storage"]  = storage
        file["metadata"] = metadata
    end
end

write_metadata(storage::NoFile, metadata) = nothing

function read_metadata(path::String)
    jldopen(joinpath(path, "METADATA.jld2"), mode[r]..., IOStream) do file
        return file["storage"], file["metadata"]
    end
end

function Base.setindex!(storage::SingleFile, block, idx::Int)
    filename   = "BLOCKS.jld2"
    objectname = @sprintf("%04d", idx)
    jldopen(joinpath(storage.path, filename), mode[a]..., IOStream) do file
        file[objectname] = block
    end
    block
end

function Base.getindex(storage::SingleFile, idx::Int)
    filename   = "BLOCKS.jld2"
    objectname = @sprintf("%04d", idx)
    jldopen(joinpath(storage.path, filename), mode[r]..., IOStream) do file
        return file[objectname]
    end
end

function Base.setindex!(storage::SingleFile, block, idx::Int, jdx::Int)
    filename   = "BLOCKS.jld2"
    objectname = @sprintf("%04d/%04d", jdx, idx)
    jldopen(joinpath(storage.path, filename), mode[a]..., IOStream) do file
        file[objectname] = block
    end
    block
end

function Base.getindex(storage::SingleFile, idx::Int, jdx::Int)
    filename   = "BLOCKS.jld2"
    objectname = @sprintf("%04d/%04d", jdx, idx)
    jldopen(joinpath(storage.path, filename), mode[r]..., IOStream) do file
        return file[objectname]
    end
end

function Base.setindex!(storage::MultipleFiles, block, idx::Int)
    filename   = @sprintf("%04d.jld2", idx)
    objectname = "block"
    jldopen(joinpath(storage.path, filename), mode[w]..., IOStream) do file
        file[objectname] = block
    end
    block
end

function Base.getindex(storage::MultipleFiles, idx::Int)
    filename   = @sprintf("%04d.jld2", idx)
    objectname = "block"
    jldopen(joinpath(storage.path, filename), mode[r]..., IOStream) do file
        return file[objectname]
    end
end

function Base.setindex!(storage::MultipleFiles, block, idx::Int, jdx::Int)
    dirname    = @sprintf("%04d",      jdx)
    filename   = @sprintf("%04d.jld2", idx)
    objectname = "block"
    isdir(joinpath(storage.path, dirname)) || mkpath(joinpath(storage.path, dirname))
    jldopen(joinpath(storage.path, dirname, filename), mode[w]..., IOStream) do file
        file[objectname] = block
    end
    block
end

function Base.getindex(storage::MultipleFiles, idx::Int, jdx::Int)
    dirname    = @sprintf("%04d",      jdx)
    filename   = @sprintf("%04d.jld2", idx)
    objectname = "block"
    jldopen(joinpath(storage.path, dirname, filename), mode[r]..., IOStream) do file
        return file[objectname]
    end
end

function Base.getindex(storage::HierarchicalStorage, m, β)
    hierarchy = storage.hierarchy

    # load each hierarchical component of the matrix
    blocks = Matrix{Complex128}[]
    dirname  = @sprintf("%04d",      β)
    filename = @sprintf("%04d.jld2", m)
    jldopen(joinpath(storage.path, dirname, filename), mode[r]..., IOStream) do file
        for idx = 1:length(hierarchy.divisions)-1
            lmax = hierarchy.divisions[idx+1]
            lmax ≥ m || continue
            objectname = @sprintf("%04d", lmax)
            push!(blocks, file[objectname])
        end
    end

    # stitch the components together into a single matrix
    output = zeros(Complex128,
                   sum(    size(block, 1) for block in blocks),
                   maximum(size(block, 2) for block in blocks))
    offset = 1
    for block in blocks
        range1 = offset:offset+size(block, 1)-1
        range2 = 1:size(block, 2)
        output_view = @view output[range1, range2]
        copy!(output_view, block)
        offset += size(block, 1)
    end
    output
end

function Base.setindex!(storage::HierarchicalStorage, block, lmax, m, β)
    dirname    = @sprintf("%04d",      β)
    filename   = @sprintf("%04d.jld2", m)
    objectname = @sprintf("%04d",   lmax)
    isdir(joinpath(storage.path, dirname)) || mkpath(joinpath(storage.path, dirname))
    jldopen(joinpath(storage.path, dirname, filename), mode[a]..., IOStream) do file
        file[objectname] = block
    end
    block
end

@inline Base.getindex(storage::Mechanism, idx::Tuple) = storage[idx...]
@inline Base.setindex!(storage::Mechanism, block, idx::Tuple) = storage[idx...] = block

function rm_old_blocks!(storage::Mechanism)
    # check to see if this path at least has a METADATA.jld2 files to help prevent this function
    # from being used to delete root or home or something catastrophic like that
    isfile(joinpath(storage.path, "METADATA.jld2")) || return
    # now delete all of the JLD2 files
    for (root, directories, files) in walkdir(storage.path)
        for file in files
            if endswith(file, ".jld2")
                rm(joinpath(root, file))
            end
        end
    end
    # finally delete all of the empty directories
    for (root, directories, files) in walkdir(storage.path, topdown=false)
        if length(files) == 0
            rm(root)
        end
    end
end

struct Cache{T}
    used  :: Ref{Bool}
    cache :: Vector{T}
end

function Cache{T}() where T
    Cache{T}(Ref(false), Vector{T}(0))
end

Base.length(cache::Cache) = length(cache.cache)
Base.getindex(cache::Cache, idx) = cache.cache[idx]
Base.setindex!(cache::Cache, X, idx) = cache.cache[idx] = copy(X)
used(cache::Cache) = cache.used[]
set!(cache::Cache) = (cache.used[] = true; cache)
unset!(cache::Cache) = (cache.used[] = false; cache)
resize!(cache::Cache, length) = Base.resize!(cache.cache, length)

