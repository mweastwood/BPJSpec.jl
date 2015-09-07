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
const HI = 1420.40575177 # MHz
const Jy = 1e-26 # one Jansky in mks unints

# cosmology
import Cosmology
const cosmology = Cosmology.cosmology()
comoving_distance(z) = Cosmology.comoving_radial_dist_mpc(cosmology,z)
age(z) = Cosmology.age_gyr(cosmology,z)
frequency(z) = HI/(1+z)
redshift(ν)  = HI/ν-1

