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

@doc """
This is a container that stores some basic information about
the interferometer and thed data being processed.
""" ->
immutable ObsParam
    freq::Vector{Float64}
    Nfreq::Int
    Nant::Int
    Nbase::Int
    Ntime::Int
end

function ObsParam(ms::Table,mmax)
    # Populate the observational parameters from the information
    # attached to a measurement set and the desired mmax to use
    # in m-mode analysis.
    lock(ms)
    freq  = MeasurementSets.frequency(ms)
    Nfreq = length(freq)
    Nant  = numrows(Table(ms[kw"ANTENNA"]))
    Nbase = div(Nant*(Nant-1),2) # no autocorrelations!
    Ntime = 2mmax+1
    unlock(ms)
    ObsParam(freq,Nfreq,Nant,Nbase,Ntime)
end

freq(obs::ObsParam) = obs.freq
Nfreq(obs::ObsParam) = obs.Nfreq
Nant(obs::ObsParam) = obs.Nant
Nbase(obs::ObsParam) = obs.Nbase
Ntime(obs::ObsParam) = obs.Ntime

@doc """
The data coming off the correlator is grouped by time.
We would instead like to have the data grouped by frequency
so that we can Fourier transform each visibility with
respect to time. This requires a corner turn.
""" ->
function cornerturn(mslist::Vector{Table},
                    obs::ObsParam;
                    outputdir::ASCIIString = "")
    @time data_by_channel, weights_by_channel = touch(obs,outputdir)
    @time for ms in mslist
        grid!(data_by_channel,weights_by_channel,ms,obs,outputdir)
    end
    @time for β = 1:Nfreq(obs)
        file = joinpath(outputdir,filename(freq(obs)[β]))
        write_data(file,data_by_channel[β],weights_by_channel[β])
    end
end

@doc """
If the output files don't exist, initialize them to zero.
Returns the current set of gridded data and weights
(ie. zero if the files don't already exist).
""" ->
function touch(obs::ObsParam,outputdir)
    data_by_channel = Matrix{Complex64}[]
    weights_by_channel = Vector{Float64}[]
    for β = 1:Nfreq(obs)
        file = joinpath(outputdir,filename(freq(obs)[β]))
        if isfile(file)
            data, weights = read_data(file)
        else
            data = zeros(Complex64,Nbase(obs),Ntime(obs))
            weights = zeros(Float64,Ntime(obs))
            write_data(file,data,weights)
        end
        push!(data_by_channel,data)
        push!(weights_by_channel,weights)
    end
    data_by_channel,weights_by_channel
end

@doc """
Grid the data in LST.
""" ->
function grid!(data_by_channel, weights_by_channel,
               ms::Table,obs::ObsParam,outputdir)
    lock(ms)

    pos   = MeasurementSets.position(ms)
    time  = MeasurementSets.time(ms)
    frame = ReferenceFrame()
    set!(frame,pos)
    lst = days(measure(frame,time,Measures.LAST))
    lst = mod(lst,1) # Discard the information on which day it is

    data  = MeasurementSets.corrected_data(ms)
    flags = MeasurementSets.flags(ms)

    grid!(data_by_channel, weights_by_channel,
          lst,data,flags,obs,outputdir)

    unlock(ms)
end

function grid!(data_by_channel, weights_by_channel,
               lst,data,flags,obs,outputdir)
    # Identify the two nearest time gridpoints
    # and calculate their corresponding weights
    pseudo_idx = lst*Ntime(obs)+1
    idx1 = floor(Int,pseudo_idx)
    idx2 =  ceil(Int,pseudo_idx)
    weight1 =  (pseudo_idx-idx1)/(idx2-idx1)
    weight2 = -(pseudo_idx-idx2)/(idx2-idx1)
    idx2 = mod(idx2-1,Ntime(obs))+1 # account for wrap-around with mod

    for β = 1:Nfreq(obs)
        gridded_data = data_by_channel[β]
        weights = weights_by_channel[β]
        # TODO: verify that Ntime and Nbase are consistent with
        # the sizes of gridded_data and weights
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

filename(ν) = @sprintf("%.3fMHz.pjd",ν/1e6)

function write_data(filename,data,weights)
    Nbase,Ntime = size(data)
    open(filename,"w") do f
        write(f,Nbase)
        write(f,Ntime)
        write(f,data)
        write(f,weights)
    end
end

function read_data(filename)
    local data, weights
    open(filename,"r") do f
        Nbase   = read(f,Int)
        Ntime   = read(f,Int)
        data    = read(f,Complex64,(Nbase,Ntime))
        weights = read(f,Float64,Ntime)
    end
    data, weights
end

