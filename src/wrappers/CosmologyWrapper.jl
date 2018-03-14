# Copyright (c) 2015-2017 Michael Eastwood
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

module CosmologyWrapper

export comoving_distance, age, frequency, redshift, approximate

using Cosmology
using Unitful, UnitfulAstro

const HI = 1420.40575177u"MHz"
const cosmology = Cosmology.FlatLCDM{Float64}(0.69,                # h
                                              0.7099122024007928,  # ΩΛ
                                              0.29,                # Ωm
                                              8.77975992071536e-5) # Ωr

"""
    comoving_distance(z)

Calculate the comoving distance (in units of Mpc) to the redshift `z`.

**Usage:**

```jldoctest
julia> comoving_distance(1)
3371.509961954628 Mpc

julia> comoving_distance(10)
9689.514711746533 Mpc
```
"""
comoving_distance(z) = Cosmology.comoving_radial_dist_mpc(cosmology, z) * u"Mpc"
comoving_distance(ν::Unitful.Frequency) = comoving_distance(redshift(ν))

function approximate(::typeof(comoving_distance), zmin, zmax)
    f = Fun(z -> Cosmology.comoving_radial_dist_mpc(cosmology, z), zmin..zmax)
    z -> f(z)*u"Mpc"
end

"""
    age(z)

Calculate the age of the universe (in units of Gyr) to the redshift `z`.

**Usage:**

```jldoctest
julia> age(1)
5.918077173774152 Gyr

julia> age(10)
0.4785005773464139 Gyr
```
"""
age(z) = Cosmology.age_gyr(cosmology, z) * u"Gyr"

"""
    frequency(z)

Calculate the frequency of the 21 cm line of Hydrogen at the redshift `z`.

**Usage:**

```jldoctest
julia> frequency(1)
710.202875885 MHz

julia> frequency(10)
129.12779561545454 MHz
```
"""
frequency(z) = HI/(1+z)

"""
    redshift(ν)

Calculate the redshift from which the emission originates if the 21 cm line
is observed at the frequency `ν`.

**Usage:**

```jldoctest
julia> redshift(100u"MHz")
13.2040575177

julia> redshift(200u"MHz")
6.10202875885
```
"""
redshift(ν)  = HI/ν-1

end

