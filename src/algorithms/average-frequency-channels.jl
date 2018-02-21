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

function average_frequency_channels(transfermatrix, Navg)
    path = dirname(transfermatrix.path)
    Nfreq  = length(transfermatrix.frequencies)
    Nfreq′ = Nfreq ÷ Navg + 1
    mmax   = transfermatrix.mmax

    # Decide which channels to average down and the resulting average frequencies.
    ranges = Array{UnitRange{Int}}(Nfreq′)
    frequencies = Array{eltype(transfermatrix.frequencies)}(Nfreq′)
    bandwidth   = Array{eltype(transfermatrix.bandwidth)}(Nfreq′)
    for idx = 1:Nfreq′
        range = (1:Navg) + (idx-1)*Navg
        range = range[1]:min(range[end], Nfreq)
        ranges[idx] = range
        frequencies[idx] = mean(transfermatrix.frequencies[range])
        bandwidth[idx]   =  sum(transfermatrix.bandwidth[range])
    end

    suffix = "-averaged"
    file = joinpath(path, "transfer-matrix"*suffix)
    output_transfermatrix = TransferMatrix(file, mmax, frequencies, bandwidth,
                                           progressbar=true, distribute=false)

    # Perform the averaging.
    prg = Progress(Nfreq′*(mmax+1))
    for idx = 1:Nfreq′
        range = ranges[idx]
        for m = 0:mmax
            B = transfermatrix[m, range[1]]
            for β in range[2:end]
                B .+= transfermatrix[m, β]
            end
            output_transfermatrix[m, idx] = B
            next!(prg)
        end
    end

    output_transfermatrix
end

