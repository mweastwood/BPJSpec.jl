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

"""
    cornerturn(mslist::Vector{Table},obs::ObsParam)

The data coming off the correlator is grouped by time.
We would instead like to have the data grouped by frequency
so that we can Fourier transform each visibility with
respect to time. This requires a corner turn.
"""
function cornerturn(mslist::Vector{Table},obs::ObsParam)
    Nms = length(mslist)
    data    = Array{Matrix{Complex64}}(Nfreq(obs))
    weights = Array{Vector{Float64}}(Nfreq(obs))
    for β = 1:Nfreq(obs)
        data[β]    = zeros(Complex64,Nbase(obs),Ntime(obs))
        weights[β] = zeros(Float64,Ntime(obs))
    end
    for ms in mslist
        grid!(data,weights,ms,obs)
    end
    for β = 1:Nfreq(obs)
        normalize!(data[β],weights[β])
    end
    data
end

"""
Grid the data in sidereal time.
"""
function grid!(data_by_channel,
               weights_by_channel,
               ms::Table,
               obs::ObsParam)
    lock(ms)

    pos   = MeasurementSets.position(ms)
    time  = MeasurementSets.time(ms)
    frame = ReferenceFrame()
    set!(frame,pos)
    set!(frame,time)
    zenith = Direction(Measures.AZEL,Quantity(0.0,Degree),Quantity(90.0,Degree))
    zenith_app  = measure(frame,zenith,Measures.APP)
    zenith_itrf = measure(frame,zenith,Measures.ITRFDIR)
    time = longitude(zenith_app,Quanta.Degree) - longitude(zenith_itrf,Quanta.Degree)
    time = mod(time/360,1)

    data  = MeasurementSets.corrected_data(ms)
    flags = MeasurementSets.flags(ms)

    grid!(data_by_channel,
          weights_by_channel,
          time,data,flags,obs)

    unlock(ms)
end

function grid!(data_by_channel,
               weights_by_channel,
               time,data,flags,obs)
    # Identify the two nearest time gridpoints
    # and calculate their corresponding weights
    pseudo_idx = time*Ntime(obs)+1
    idx1 = floor(Int,pseudo_idx)
    idx2 =  ceil(Int,pseudo_idx)
    weight1 =  (pseudo_idx-idx1)/(idx2-idx1)
    weight2 = -(pseudo_idx-idx2)/(idx2-idx1)
    idx2 = mod(idx2-1,Ntime(obs))+1 # account for wrap-around with mod

    for β = 1:Nfreq(obs)
        gridded_data = data_by_channel[β]
        weights = weights_by_channel[β]
        weights[idx1] += weight1
        weights[idx2] += weight2
        α1 = 1 # baseline index into data (counts autocorrelations)
        α2 = 1 # baseline index into gridded_data (drops autocorrelations)
        for ant1 = 1:Nant(obs), ant2 = ant1:Nant(obs)
            if ant1 != ant2
                if !flags[1,β,α1] && !flags[4,β,α1]
                    gridded_data[α2,idx1] += weight1*(data[1,β,α1]+data[4,β,α1])/2
                    gridded_data[α2,idx2] += weight2*(data[1,β,α1]+data[4,β,α1])/2
                end
                α2 += 1
            end
            α1 += 1
        end
    end
end

function normalize!(data,weights)
    for idx = 1:size(data,2), α = 1:size(data,1)
        data[α,idx] /= weights[idx]
    end
end

