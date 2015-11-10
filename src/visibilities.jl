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
    create_empty_visibilities(filename, Nbase, Ntime, frequencies)

Create the JLD file that will store the data that
feeds into $m$-mode analysis.
"""
function create_empty_visibilities(filename, Nbase, Ntime, frequencies)
    Nfreq = length(frequencies)
    jldopen(filename,"w",compress=true) do file
        file["Nfreq"] = Nfreq
        file["Nbase"] = Nbase
        file["Ntime"] = Ntime
        for β = 1:Nfreq
            name  = @sprintf("%.3f",frequencies[β]/1e6)
            group = g_create(file,name)
            group["data"] = zeros(Complex128,Nbase,Ntime)
            group["weights"] = zeros(Float64,Nbase,Ntime)
        end
    end
end

"""
    grid_visibilities(filename, ms::MeasurementSet)

Grid the data from the measurement set.
"""
function grid_visibilities(filename, ms::MeasurementSet)
    isfile(filename) || error("$filename does not exist!")

    # Note that the following value for Nbase does not count the autocorrelations.
    # On the other hand ms.Nbase does count the autocorrelations. These two
    # quantities are not equal.
    Nbase = div(ms.Nant*(ms.Nant-1),2)

    time  = sidereal_time(ms)
    data  = TTCal.get_corrected_data(ms)
    flags = TTCal.get_flags(ms)

    # keep only the "Stokes I" visibilities and discard the
    # autocorrelations because they have an additive noise bias
    reduced_data  = zeros(Complex64,Nbase,ms.Nfreq)
    reduced_flags = zeros(     Bool,Nbase,ms.Nfreq)

    idx = 1
    for α = 1:ms.Nbase
        ms.ant1[α] == ms.ant2[α] && continue
        for β = 1:ms.Nfreq
            reduced_data[idx,β]  = 0.5*(data[1,β,α] + data[4,β,α])
            reduced_flags[idx,β] = flags[1,β,α] || flags[4,β,α]
        end
        idx += 1
    end

    grid_visibilities(filename, reduced_data, reduced_flags, frequencies, time)
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
    mod(time/360,1)
end

doc"""
    grid_visibilities(filename, data, flags, frequencies, time)

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
function grid_visibilities(filename, data, flags, frequencies, time)
    jldopen(filename,"r+",compress=true) do file
        Nfreq = file["Nfreq"] |> read
        Nbase = file["Nbase"] |> read
        Ntime = file["Ntime"] |> read

        # Identify the two nearest time gridpoints
        # and calculate their corresponding weights
        Δt = 1/Ntime
        grid = 0.0:Δt:(1.0-Δt)
        idx1 = searchsortedlast(grid,time)
        idx2 = idx1 == Ntime? 1 : idx1+1
        weight1 = 1 - abs(time - grid[idx1])/Δt
        weight2 = 1 - abs(time - grid[idx2])/Δt

        for β = 1:Nfreq
            name  = @sprintf("%.3f",frequencies[β]/1e6)
            group = file[name]

            # use type assertions here to prevent type uncertainty
            # from slowing down the inner loop
            data1 = group["data"][:,idx1]::Array{Complex128,2}
            data2 = group["data"][:,idx2]::Array{Complex128,2}
            weights1 = group["weights"][:,idx1]::Array{Float64,2}
            weights2 = group["weights"][:,idx2]::Array{Float64,2}

            for α = 1:Nbase
                flags[α,β] && continue
                data1[α] += weight1 * data[α,β]
                data2[α] += weight2 * data[α,β]
                weights1[α] += weight1
                weights2[α] += weight2
            end

            group["data"][:,idx1] = data1
            group["data"][:,idx2] = data2
            group["weights"][:,idx1] = weights1
            group["weights"][:,idx2] = weights2
        end
    end
end

function load_visibilities(filename, frequency)
    local data, weights
    jldopen(filename,"r") do file
        Nfreq = file["Nfreq"] |> read
        Nbase = file["Nbase"] |> read
        Ntime = file["Ntime"] |> read

        name  = @sprintf("%.3f",frequency/1e6)
        group = file[name]
        rawdata = group["data"] |> read
        weights = group["weights"] |> read
        data = reinterpret(Complex128,rawdata,(Nbase,Ntime))
    end
    for i in eachindex(data,weights)
        weights[i] == 0 && continue
        data[i] = data[i] / weights[i]
    end
    data
end

