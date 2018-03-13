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

function permute_mmodes(input, transfermatrix; storage=NoFile(), progress=false)
    output = similar(input, storage)
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
            remotecall_fetch(_permute_mmodes, pool,
                             input, output, transfermatrix, m, β)
            progress && increment()
        end
    end
    output
end

function _permute_mmodes(input, output, transfermatrix, m, β)
    permutation = baseline_permutation(transfermatrix, m)
    output[m, β] = input[m, β][permutation]
end

