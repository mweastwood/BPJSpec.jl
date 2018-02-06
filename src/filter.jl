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

function filter(transfermatrix, signal, foregrounds, noise)
    # TODO check that ' is doing the conjugate transpose and not the regular transpose
    # (@. might be rewriting ' as .' which is the regular transpose)
    @. S = B*signal*B'
    @. F = B*foreground*B'
    # compute eig(S, F) -> R filter
    # R'*S and R'*N

end








