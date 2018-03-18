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

function average_frequency_channels(input::Union{MFBlockVector, MFBlockMatrix},
                                    Navg; storage=NoFile(), progress=false)
    ν  = input.frequencies
    Δν = input.bandwidth
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
    output = similar(input, storage, input.mmax, ν′, Δν′)
    queue  = collect(indices(output))
    pool   = CachingPool(workers())
    if progress
        lck = ReentrantLock()
        prg = Progress(length(queue))
        increment() = (lock(lck); next!(prg); unlock(lck))
    end
    @sync for worker in workers()
        @async while length(queue) > 0
            m, β = shift!(queue)
            remotecall_fetch(_average_frequency_channels, pool,
                             input, output, Δν, Δν′, m, β, partition[β])
            progress && increment()
        end
    end
    output
end

function _average_frequency_channels(input, output, Δν, Δν′, m, β, channels)
    β′ = channels[1]
    weight = u(NoUnits, Δν[β′]/Δν′[β])
    B = input[m, β′] .* weight

    for β′ in channels[2:end]
        weight = u(NoUnits, Δν[β′]/Δν′[β])
        B .+= input[m, β′] .* weight
    end

    output[m, β] = B
end

