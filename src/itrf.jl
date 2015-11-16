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
    itrf_baselines(ms::MeasurementSet)

Calculate the $(u,v,w)$ coordinates in the ITRF coordinate system with the phase
center at true north.
"""
function itrf_baselines(ms::MeasurementSet)
    antenna_table = Table(ms.table[kw"ANTENNA"])
    positions = antenna_table["POSITION"]
    unlock(antenna_table)

    # calculate the baselines while skipping the autocorrelations
    Nbase = ms.Nbase - ms.Nant
    u = zeros(Nbase)
    v = zeros(Nbase)
    w = zeros(Nbase)
    idx = 1
    for α = 1:ms.Nbase
        ms.ant1[α] == ms.ant2[α] && continue
        u[idx] = positions[1,ms.ant1[α]] - positions[1,ms.ant2[α]]
        v[idx] = positions[2,ms.ant1[α]] - positions[2,ms.ant2[α]]
        w[idx] = positions[3,ms.ant1[α]] - positions[3,ms.ant2[α]]
        idx += 1
    end
    u,v,w
end

doc"""
    itrf_phasecenter(ms::MeasurementSet)

Calculate the $(l,m,n)$ coordinates corresponding to the actual phase
center of the array with respsect to true north.
"""
function itrf_phasecenter(ms::MeasurementSet)
    itrf = measure(ms.frame,ms.phase_direction,dir"ITRF")
    vector(itrf)
end

doc"""
    itrf_beam(frame::ReferenceFrame, beam::BeamModel, frequency)

Create a Healpix map of a TTCal beam model.
At the moment only the $I \rightarrow I$ element of the Mueller
matrix is used.
"""
function itrf_beam(frame::ReferenceFrame, beam::TTCal.BeamModel, frequency)
    zenith = Direction(dir"AZEL",Quantity(0.0,"deg"),Quantity(90.0,"deg"))
    zenith_itrf = measure(frame,zenith,dir"ITRF")
    zenith_vec  = [vector(zenith_itrf)...]
    global_north = [0.0,0.0,1.0]
    local_north  = gramschmidt(global_north,zenith_vec)
    map = HealpixMap(zeros(nside2npix(512)))
    @showprogress 1 "Creating map of the beam..." for i = 1:length(map)
        vec = LibHealpix.pix2vec_ring(512,i)
        el  = π/2 - angle_between(vec,zenith_vec)
        vec = gramschmidt(vec,zenith_vec)
        az  = angle_between(vec,local_north)
        if el < 0
            map[i] = 0
        else
            J = beam(frequency,az,el)
            M = MuellerMatrix(J)
            map[i] = M.mat[1,1]
        end
    end
    map
end

