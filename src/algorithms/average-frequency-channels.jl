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

function average_frequency_channels(transfermatrix, Navg; output="")
    path = dirname(transfermatrix.path)
    Nfreq  = length(transfermatrix.frequencies)
    Nfreq′ = cld(Nfreq, Navg)
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

    if output == ""
        suffix = "-averaged"
        file = joinpath(path, "transfer-matrix"*suffix)
    else
        file = output
    end
    output_transfermatrix = TransferMatrix(file, mmax, frequencies, bandwidth,
                                           progressbar=true, distribute=true)

    # Perform the averaging.
    queue = [(m, β) for β = 1:Nfreq′ for m = 0:mmax]
    pool  = CachingPool(workers())
    lck = ReentrantLock()
    prg = Progress(length(queue))
    increment() = (lock(lck); next!(prg); unlock(lck))
    @sync for worker in workers()
        @async while length(queue) > 0
            m, β = shift!(queue)
            remotecall_fetch(_average_frequency_channels, pool,
                             transfermatrix, output_transfermatrix, m, β, ranges[β])
            increment()
        end
    end

    output_transfermatrix
end

function _average_frequency_channels(input, output, m, β, range)
    β′ = range[1]
    weight = uconvert(NoUnits, input.bandwidth[β′]/output.bandwidth[β])
    B = input[m, β′] .* weight

    for β′ in range[2:end]
        weight = uconvert(NoUnits, input.bandwidth[β′]/output.bandwidth[β])
        B .+= input[m, β′] .* weight
    end

    output[m, β] = B
end

