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
    freq::Vector{Float64}          # list of frequency channels
    antpos::Matrix{Float64}        # 3 by Nant array containing the (x,y,z) ITRF positions of each antenna
    phasecenter::NTuple{3,Float64} # unit vector specifying the ITRF phase center
    Nfreq::Int
    Nant::Int
    Nbase::Int
    Ntime::Int
end

function ObsParam(ms::Table;mmax::Int=100)
    # Populate the observational parameters from the information
    # attached to a measurement set and the desired mmax to use
    # in m-mode analysis.
    lock(ms)
    Ntime = 2mmax+1

    # Frequency channels
    freq  = MeasurementSets.frequency(ms)
    Nfreq = length(freq)

    # Antenna positions
    antenna_table = Table(ms[kw"ANTENNA"])
    antpos = antenna_table["POSITION"]
    Nant = size(antpos,2)
    Nbase = div(Nant*(Nant-1),2) # no autocorrelations!

    # Phase center
    frame = ReferenceFrame()
    set!(frame,MeasurementSets.time(ms))
    set!(frame,Measures.from_xyz_in_meters(Measures.ITRF,antpos[1,1],antpos[2,1],antpos[3,1]))
    phasecenter_J2000 = MeasurementSets.phase_direction(ms)
    phasecenter_ITRF  = measure(frame,phasecenter_J2000,Measures.ITRFDIR)
    phasecenter = Measures.xyz_in_meters(phasecenter_ITRF)

    unlock(ms)
    ObsParam(freq,antpos,phasecenter,Nfreq,Nant,Nbase,Ntime)
end

freq(obs::ObsParam) = obs.freq
antpos(obs::ObsParam) = obs.antpos
phasecenter(obs::ObsParam) = obs.phasecenter
Nfreq(obs::ObsParam) = obs.Nfreq
Nant(obs::ObsParam) = obs.Nant
Nbase(obs::ObsParam) = obs.Nbase
Ntime(obs::ObsParam) = obs.Ntime

