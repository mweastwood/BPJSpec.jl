# Copyright (c) 2015 Michael Eastwood
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

macro write_metadata(file,description,version)
    quote
        attrs($file)["description"] = $description
        attrs($file)["version_major"] = $version.major
        attrs($file)["version_minor"] = $version.minor
    end
end

macro check_version(file,desired_version)
    quote
        version = VersionNumber(read(attrs($file)["version_major"]),
                                read(attrs($file)["version_minor"]))
        if version != $desired_version
            error("Incorrect version $version.")
        end
    end
end

@doc """
Reinterpret a complex array as a float array.
This is necessary because HDF5 does not support
complex numbers.
""" ->
function to_float{T<:FloatingPoint}(::Type{T},x)
    sz = size(x)
    sz = (2sz[1],sz[2:end]...)
    reinterpret(T,x,sz)
end

@doc """
Reinterpret a float array as a complex array.
This is necessary because HDF5 does not support
complex numbers.
""" ->
function to_complex{T<:Complex}(::Type{T},x)
    sz = size(x)
    sz = (div(sz[1],2),sz[2:end]...)
    reinterpret(T,x,sz)
end

const DATA_VER = v"1.0"

@doc """
Write gridded visibility data to disk.
""" ->
function write_data(filename,
                    data::Matrix{Complex64},
                    weights::Vector{Float64})
    Nbase,Ntime = size(data)
    h5open(filename,"w") do file
        @write_metadata file "BPJSpec data file" DATA_VER
        attrs(file)["Nbase"] = Nbase
        attrs(file)["Ntime"] = Ntime
        file["data"] = to_float(Float32,data)
        file["weights"] = weights
    end
    nothing
end

@doc """
Read gridded visibility data from disk.
""" ->
function read_data(filename)
    local Nbase, Ntime
    local data, weights
    h5open(filename,"r") do file
        @check_version file DATA_VER
        Nbase = read(attrs(file)["Nbase"])
        Ntime = read(attrs(file)["Ntime"])
        data = to_complex(Complex64,read(file["data"]))
        weights = read(file["weights"])
    end
    data, weights
end

const MMODE_VER = v"1.0"

@doc """
Write m-modes to disk.
""" ->
function write_mmodes(filename,
                      v::MModes)
    h5open(filename,"w") do file
        @write_metadata file "BPJSpec m-modes" MMODE_VER
        attrs(file)["Nbase"] = Nbase(v)
        attrs(file)["mmax"] = mmax(v)
        g_create(file,"blocks")
        blocks_group = file["blocks"]
        for m = 0:mmax(B)
            blocks_group[string(m)] = to_float(Float64,block(v,m))
        end
    end
    nothing
end

@doc """
Read m-modes from disk.
""" ->
function read_mmodes(filename)
    local Nbase, mmax
    local blocks
    h5open(filename,"r") do file
        @check_version file MMODE_VER
        Nbase = read(attrs(file)["Nbase"])
        mmax = read(attrs(file)["mmax"])
        blocks_group = file["blocks"]
        blocks = Matrix{Complex128}[to_complex(Complex128,read(blocks_group[string(m)]))
                                        for m = 0:mmax]
    end
    MModes{Nbase,mmax}(blocks)
end

const TRANSFER_VER = v"1.0"

@doc """
Write a transfer matrix to disk.
""" ->
function write_transfermatrix(filename,
                              B::TransferMatrix)
    h5open(filename,"w") do file
        @write_metadata file "BPJSpec transfer matrix" TRANSFER_VER
        attrs(file)["Nbase"] = Nbase(B)
        attrs(file)["lmax"] = lmax(B)
        attrs(file)["mmax"] = mmax(B)
        g_create(file,"blocks")
        blocks_group = file["blocks"]
        for m = 0:mmax(B)
            blocks_group[string(m)] = to_float(Float64,block(B,m))
        end
    end
    nothing
end

@doc """
Read a transfer matrix from disk.
""" ->
function read_transfermatrix(filename)
    local description, version
    local Nbase, lmax, mmax
    local blocks
    h5open(filename,"r") do file
        @check_version file TRANSFER_VER
        Nbase = read(attrs(file)["Nbase"])
        lmax = read(attrs(file)["lmax"])
        mmax = read(attrs(file)["mmax"])
        blocks_group = file["blocks"]
        blocks = Matrix{Complex128}[to_complex(Complex128,read(blocks_group[string(m)]))
                                        for m = 0:mmax]
    end
    TransferMatrix{Nbase,lmax,mmax}(blocks)
end

