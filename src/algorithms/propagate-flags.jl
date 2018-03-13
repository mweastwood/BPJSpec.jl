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

function propagate_flags(input_mmodes, input_transfermatrix;
                         mmodes_storage=NoFile(),
                         transfermatrix_storage=NoFile(),
                         progress=false)

    output_mmodes         = similar(input_mmodes,                 mmodes_storage)
    output_transfermatrix = similar(input_transfermatrix, transfermatrix_storage)
    multi_broadcast!(_propagate_flags,
                     (output_mmodes, output_transfermatrix),
                     (input_mmodes, input_transfermatrix), progress=progress)
end

function _propagate_flags(mmodes, transfermatrix)
    notflagged = mmodes .!= 0
    mmodes[notflagged], transfermatrix[notflagged, :]
end

