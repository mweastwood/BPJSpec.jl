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

# constants
const c = 2.99792e8 # m/s
const k = 1.38065e-23 # J/K
const HI = 1420.40575177e6 # Hz
const Jy = 1e-26 # one Jansky in mks unints

# cosmology
const cosmology = Cosmology.cosmology()

"""
    comoving_distance(z)

Calculate the comoving distance (in units of Mpc) to the redshift `z`.
"""
comoving_distance(z) = Cosmology.comoving_radial_dist_mpc(cosmology,z)

"""
    age(z)

Calculate the age of the universe (in units of Gyr) to the redshift `z`.
"""
age(z) = Cosmology.age_gyr(cosmology,z)

"""
    frequency(z)

Calculate the frequency (in Hz) of the 21 cm line of Hydrogen at the redshift `z`.
"""
frequency(z) = HI/(1+z)

"""
    redshift(ν)

Calculate the redshift from which the emission originates if the 21 cm line
is observed at the frequency `ν` (in Hz).
"""
redshift(ν)  = HI/ν-1

"""
    beam_solid_angle(beam::SineBeam)

Calculate the solid angle (in steradians) of the given beam model.
"""
function beam_solid_angle(beam::SineBeam)
    2π*quadgk(x->sin(x)^beam.α*cos(x),0,π/2)[1]
end

function jansky_to_kelvin(x, ν, Ω)
    out * (c^2/(2*ν^2*k*Ω) * Jy)
end

