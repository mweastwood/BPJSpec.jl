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

module BPJSpec

export CGTable, multiply

using HEALPix
using CasaCore.Quanta
using CasaCore.Measures
using CasaCore.Tables
using TTCal
using JSON

include("special.jl") # special functions
include("clebschgordan.jl")

function readbeam(filename::AbstractString)
    dict = JSON.parsefile(filename)
    alm = Alm(Complex128,dict["lmax"],dict["mmax"])
    real_part = dict["alm_real"]
    imag_part = dict["alm_imag"]
    for i = 1:length(alm)
        alm[i] = complex(real_part[i],imag_part[i])
    end
    alm
end

function writebeam(filename::AbstractString,beam::Alm)
    dict = Dict{UTF8String,Any}()
    dict["lmax"] = lmax(beam)
    dict["mmax"] = mmax(beam)
    alm = coefficients(beam)
    dict["alm_real"] = real(alm)
    dict["alm_imag"] = imag(alm)
    file = open(filename,"w")
    JSON.print(file,dict)
    close(file)
    nothing
end

#include("tests.jl")

end

