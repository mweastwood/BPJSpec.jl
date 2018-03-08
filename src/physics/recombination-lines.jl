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

const R∞  = 1.0973731568508e7u"m^-1" # Rydberg constant
const amu = 1.66053886e-27u"kg"      # atomic mass unit

function radio_recombination_line(nuclear_mass, n, Δn)
    # NOTE: this function is adapted from a Matlab routine provided by Judd Bowman.
    Z = 1 # the nuclear charge, for RRLs we can almost always assume Z=1 since atoms can be
          # considered hydrogenic in that the nuclear charge is screened by electrons in the low
          # energy states
    Rx = R∞ / (1 + u"me"/nuclear_mass)
    uconvert(u"MHz", Z^2 * Rx * u"c" * (1/n^2 - 1/(n+Δn)^2))
end

hydrogen(n, Δn) = radio_recombination_line(1.007825amu, n, Δn)
carbon(n, Δn)   = radio_recombination_line(      12amu, n, Δn)

