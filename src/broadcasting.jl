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

progressbar(matrix) = matrix.progressbar
distribute(matrix) = matrix.distribute
#flatten(x) = reshape(x, length(x))

function Base.broadcast!(f, output::Union{BlockVector, BlockMatrix}, args...)
    if distribute(output)
        return distributed_broadcast!(f, output, args...)
    else
        return local_broadcast!(f, output, args...)
    end
end

function local_broadcast!(f, output, args...)
    #queue = flatten(collect(Iterators.product(indices(output)...)))
    queue = indices(output)
    if progresbar(output)
        prg = Progress(length(queue))
    end
    for indices in queue
        just_do_it!(f, output, args, indices)
    end
    output
end

function distributed_broadcast!(f, output, args...)
    #queue = flatten(collect(Iterators.product(indices(output)...)))
    queue = indices(output)
    pool  = CachingPool(workers())
    if progressbar(output)
        lck = ReentrantLock()
        prg = Progress(length(queue))
        increment() = (lock(lck); next!(prg); unlock(lck))
    end
    @sync for worker in workers()
        @async while length(queue) > 0
            indices = shift!(queue)
            remotecall_fetch(just_do_it!, pool, f, output, args, indices)
            progressbar(output) && increment()
        end
    end
    output
end

function just_do_it!(f, output, args, indices) # (c) Nike
    output[indices...] = f(getindex.(args, indices...)...)
end

