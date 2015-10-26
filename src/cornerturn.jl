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

doc"""
    touch(filename, Nfreq, Nbase, Ntime)

Create the HDF5 file that will store the data that
feeds into $m$-mode analysis.
"""
function touch(filename, Nfreq, Nbase, Ntime)
    h5open(filename,"w") do file
        attrs(file)["Nfreq"] = Nfreq
        attrs(file)["Nbase"] = Nbase
        attrs(file)["Ntime"] = Ntime
        for β = 1:Nfreq
            group = g_create(file,string(β))
            # create datasets to store the real and imaginary
            # components of the visibilities
            group["real"] = zeros(Nbase,Ntime)
            group["imag"] = zeros(Nbase,Ntime)
            # create a dataset that essentially counts the fractional
            # number of integrations that go into each sidereal time
            group["weights"] = zeros(Nbase,Ntime)
        end
    end
end

function grid(filename, ms::MeasurementSet)
    isfile(filename) || error("$filename does not exist!")
    time  = sidereal_time(ms)
    data  = TTCal.get_corrected_data(ms)
    flags = TTCal.get_flags(ms)

    # keep only the "Stokes I" visibilities and discard the
    # autocorrelations because they have an additive noise bias
    reduced_Nbase = div(ms.Nant*(ms.Nant-1),2)
    reduced_data  = zeros(Complex64,reduced_Nbase,ms.Nfreq)
    reduced_flags = zeros(     Bool,reduced_Nbase,ms.Nfreq)

    idx = 1
    for α = 1:ms.Nbase
        ms.ant1[α] == ms.ant2[α] && continue
        for β = 1:ms.Nbase
            reduced_data[idx,β]  = 0.5*(data[1,β,α] + data[4,β,α])
            reduced_flags[idx,β] = flags[1,β,α] || flags[4,β,α]
        end
        idx += 1
    end

    grid(filename, time, reduced_data, reduced_flags)
end

doc"""
    sidereal_time(ms::MeasurementSet)

Get the local sidereal time on the interval $\left(0,1\right]$ for
the given measurement set.
"""
function sidereal_time(ms::MeasurementSet)
    zenith = Direction(dir"AZEL",Quantity(0.0,"deg"),Quantity(90.0,"deg"))
    zenith_app = measure(ms.frame,zenith,Measures.APP)
    time = longitude(zenith_app,Quanta.Degree)
    mod(time/360,1)
end

doc"""
    grid(filename, time, data, flags)

For a sidereal time on the interval $\left(0,1\right]$, grid and
add the data to the given file.

Here "grid" means to put the data onto a set of points that are
evenly space in sidereal time. If the correlator is continuously
running, it will spit out an integration every 30 seconds (for
example). However 30 seconds does not divide evenly into a sidereal
day so that these integrations will slowly start slipping with
respect to sidereal time as we begin to average more days together.
By gridding onto a sidereal time grid, we avoid this problem.
"""
function grid(filename, time, data, flags)
    h5open(filename,"r+") do file
        Nfreq = attrs(file)["Nfreq"] |> read
        Nbase = attrs(file)["Nbase"] |> read
        Ntime = attrs(file)["Ntime"] |> read

        # Identify the two nearest time gridpoints
        # and calculate their corresponding weights
        Δt = 1/Ntime
        grid = 0.0:Δt:(1.0-Δt)
        idx1 = searchsortedlast(grid,time)
        idx2 = idx1 == Ntime? 0 : idx1+1
        weight1 = 1 - abs(time - grid[idx1])/Δt
        weight2 = 1 - abs(time - grid[idx2])/Δt

        for β = 1:Nfreq
            group = file[string(β)]
            realpart = group["real"][:,idx1:idx2]::Matrix{Float64}
            imagpart = group["imag"][:,idx1:idx2]::Matrix{Float64}
            weights  = group["weights"][:,idx1:idx2]::Matrix{Float64}

            for α = 1:Nbase
                flags[α,β] && continue
                realpart[α,1] += weight1 * real(data[α,β])
                imagpart[α,1] += weight1 * imag(data[α,β])
                realpart[α,2] += weight2 * real(data[α,β])
                imagpart[α,2] += weight2 * imag(data[α,β])
                weights[α,1] += weight1
                weights[α,2] += weight2
            end

            group["real"][:,idx1:idx2]  = realpart
            group["imag"][:,idx1:idx2]  = imagpart
            group["weights"][:,idx1:idx2] = weights
        end
    end
end

