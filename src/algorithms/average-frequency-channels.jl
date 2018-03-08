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

function average_frequency_channels(input_matrix, Navg; output=NoFile(), progress=false)
    ν = frequencies(input_matrix)
    Δν = bandwidth(input_matrix)
    Nfreq = length(ν)

    partition = collect(Iterators.partition(1:Nfreq, Navg))
    weights   = u.(u"Hz", Δν)

    Nfreq′ = length(partition)
    ν′  = similar(ν,  Nfreq′)
    Δν′ = similar(Δν, Nfreq′)
    for β = 1:length(partition)
        channels = partition[β]
        ν′[β]  = sum(weights[channels].*ν[channels]) / sum(weights[channels])
        Δν′[β] = sum(Δν[channels])
    end

    # Perform the averaging.
    output_matrix = MFBlockMatrix(output, mmax(input_matrix), ν′, Δν′)
    queue = collect(indices(output_matrix))
    pool  = CachingPool(workers())
    if progress
        lck = ReentrantLock()
        prg = Progress(length(queue))
        increment() = (lock(lck); next!(prg); unlock(lck))
    end
    @sync for worker in workers()
        @async while length(queue) > 0
            m, β = shift!(queue)
            remotecall_fetch(_average_frequency_channels, pool,
                             input_matrix, output_matrix, m, β, partition[β])
            progress && increment()
        end
    end

    output_matrix
end

function _average_frequency_channels(input, output, m, β, channels)
    β′ = channels[1]
    weight = u(NoUnits, bandwidth(input)[β′]/bandwidth(output)[β])
    B = input[m, β′] .* weight

    for β′ in channels[2:end]
        weight = u(NoUnits, bandwidth(input)[β′]/bandwidth(output)[β])
        B .+= input[m, β′] .* weight
    end

    output[m, β] = B
end

