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

