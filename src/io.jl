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

const DATA_VER = v"1.0"

@doc """
Write gridded visibility data to disk.
""" ->
function write_data(filename,
                    data::Matrix{Complex64},
                    weights::Vector{Float64})
    Nbase,Ntime = size(data)
    h5open(filename,"w") do file
        attrs(file)["description"] = "BPJSpec data file"
        attrs(file)["version_major"] = DATA_VER.major
        attrs(file)["version_minor"] = DATA_VER.minor
        attrs(file)["Nbase"] = Nbase
        attrs(file)["Ntime"] = Ntime
        # HDF5 doesn't support complex numbers,
        # so we circumvent this by reinterpreting
        # the complex numbers as two floating point
        # numbers.
        file["data"] = reinterpret(Float32,data,(2Nbase,Ntime))
        file["weights"] = weights
    end
    nothing
end

@doc """
Read gridded visibility data from disk.
""" ->
function read_data(filename)
    local description, version
    local Nbase, Ntime
    local data, weights
    h5open(filename,"r") do file
        description = read(attrs(file)["description"])
        version = VersionNumber(read(attrs(file)["version_major"]),
                                read(attrs(file)["version_minor"]))
        if version != DATA_VER
            error("Incorrect version $version.")
        end
        Nbase = read(attrs(file)["Nbase"])
        Ntime = read(attrs(file)["Ntime"])
        data = reinterpret(Complex64,read(file["data"]),(Nbase,Ntime))
        weights = read(file["weights"])
    end
    data, weights
end

const TRANSFER_VER = v"1.0"

@doc """
Write a transfer matrix to disk.
""" ->
function write_transfermatrix(filename,
                              B::TransferMatrix)
    h5open(filename,"w") do file
        attrs(file)["description"] = "BPJSpec transfer matrix"
        attrs(file)["version_major"] = TRANSFER_VER.major
        attrs(file)["version_minor"] = TRANSFER_VER.minor
        attrs(file)["Nbase"] = Nbase(B)
        attrs(file)["lmax"] = lmax(B)
        attrs(file)["mmax"] = mmax(B)
        g_create(file,"blocks")
        blocks_group = file["blocks"]
        # HDF5 doesn't support complex numbers,
        # so we circumvent this by reinterpreting
        # the complex numbers as two floating point
        # numbers.
        for m = 0:mmax(B)
            blocks_group[string(m)] = reinterpret(Float64,block(B,m),(2*two(m)*Nbase(B),lmax(B)-m+1))
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
        description = read(attrs(file)["description"])
        version = VersionNumber(read(attrs(file)["version_major"]),
                                read(attrs(file)["version_minor"]))
        if version != DATA_VER
            error("Incorrect version $version.")
        end
        Nbase = read(attrs(file)["Nbase"])
        lmax = read(attrs(file)["lmax"])
        mmax = read(attrs(file)["mmax"])
        blocks_group = file["blocks"]
        blocks = Matrix{Complex128}[reinterpret(Complex128,
                                                read(blocks_group[string(m)]),
                                                (two(m)*Nbase,lmax-m+1)) for m = 0:mmax]
    end
    TransferMatrix{Nbase,lmax,mmax}(blocks)
end

