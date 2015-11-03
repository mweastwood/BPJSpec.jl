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

"""
    grid(filename, ms::MeasurementSet)

Grid the data from the measurement set.
"""
grid(filename, ms::MeasurementSet) = grid(filename, ms, 1:ms.Nfreq)

"""
    grid(filename, ms::MeasurementSet, channels)

Grid the data from the measurement set, but select a subset
of the frequency channels.
"""
function grid(filename, ms::MeasurementSet, channels)
    isfile(filename) || error("$filename does not exist!")
    time  = sidereal_time(ms)
    data  = TTCal.get_corrected_data(ms)
    flags = TTCal.get_flags(ms)
    Nfreq = length(channels)

    # keep only the "Stokes I" visibilities and discard the
    # autocorrelations because they have an additive noise bias
    reduced_Nbase = div(ms.Nant*(ms.Nant-1),2)
    reduced_data  = zeros(Complex64,reduced_Nbase,ms.Nfreq)
    reduced_flags = zeros(     Bool,reduced_Nbase,ms.Nfreq)

    idx = 1
    for α = 1:ms.Nbase
        ms.ant1[α] == ms.ant2[α] && continue
        for β = 1:Nfreq
            reduced_data[idx,β]  = 0.5*(data[1,channels[β],α] + data[4,channels[β],α])
            reduced_flags[idx,β] = flags[1,channels[β],α] || flags[4,channels[β],α]
        end
        idx += 1
    end

    grid(filename, time, reduced_data, reduced_flags)
end

doc"""
    sidereal_time(ms::MeasurementSet)

Get the local sidereal time on the interval $0 \le t < 1$ for
the given measurement set.
"""
function sidereal_time(ms::MeasurementSet)
    zenith = Direction(dir"AZEL",Quantity(0.0,"deg"),Quantity(90.0,"deg"))
    zenith_app = measure(ms.frame,zenith,dir"APP")
    time = longitude(zenith_app,"deg")
    # BEGIN TEMP FIX (until I regenerate transfer matrices)
    zenith_itrf = measure(ms.frame,zenith,dir"ITRF")
    time -= longitude(zenith_itrf,"deg")
    # END TEMP FIX
    mod(time/360,1)
end

doc"""
    grid(filename, time, data, flags)

For a sidereal time on the interval $0 \le t < 1$, grid and
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
        idx2 = idx1 == Ntime? 1 : idx1+1
        weight1 = 1 - abs(time - grid[idx1])/Δt
        weight2 = 1 - abs(time - grid[idx2])/Δt

        for β = 1:Nfreq
            group = file[string(β)]

            realpart1 = group["real"][:,idx1]::Matrix{Float64}
            imagpart1 = group["imag"][:,idx1]::Matrix{Float64}
            realpart2 = group["real"][:,idx2]::Matrix{Float64}
            imagpart2 = group["imag"][:,idx2]::Matrix{Float64}
            weights1  = group["weights"][:,idx1]::Matrix{Float64}
            weights2  = group["weights"][:,idx2]::Matrix{Float64}

            for α = 1:Nbase
                flags[α,β] && continue
                realpart1[α] += weight1 * real(data[α,β])
                imagpart1[α] += weight1 * imag(data[α,β])
                realpart2[α] += weight2 * real(data[α,β])
                imagpart2[α] += weight2 * imag(data[α,β])
                weights1[α] += weight1
                weights2[α] += weight2
            end

            group["real"][:,idx1]  = realpart1
            group["imag"][:,idx1]  = imagpart1
            group["real"][:,idx2]  = realpart2
            group["imag"][:,idx2]  = imagpart2
            group["weights"][:,idx1] = weights1
            group["weights"][:,idx2] = weights2
        end
    end
end

function load_visibilities(filename, channel)
    local Nfreq, Nbase
    local data, weights
    h5open(filename,"r") do file
        Nfreq = attrs(file)["Nfreq"] |> read
        Nbase = attrs(file)["Nbase"] |> read

        group = file[string(channel)]
        realpart = group["real"] |> read
        imagpart = group["imag"] |> read
        weights  = group["weights"] |> read
        data = complex(realpart,imagpart)
    end
    for i in eachindex(data,weights)
        weights[i] == 0 && continue
        data[i] = data[i] / weights[i]
    end
    # Zero any baselines that are missing data
    #for α = 1:Nbase
    #    any(slice(weights,α,:) .== 0) || continue
    #    data[α,:] = 0
    #end
    data
end

