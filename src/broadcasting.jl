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

progressbar(matrix::BlockMatrix) = matrix.progressbar
distribute(matrix::BlockMatrix)  = matrix.distribute
progressbar(::BlockVector) = false
distribute(::BlockVector)  = false

function Base.broadcast!(f, output::Union{BlockVector, BlockMatrix}, args...)
    if distribute(output)
        return distributed_broadcast!(f, output, args)
    else
        return local_broadcast!(f, output, args)
    end
end

function local_broadcast!(f, output, args)
    queue = indices(output)
    if progressbar(output)
        prg = Progress(length(queue))
    end
    for indices in queue
        just_do_it!(f, output, args, indices)
    end
    output
end

function distributed_broadcast!(f, output, args)
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

function multi_broadcast!(f, outputs, args)
    queue = indices(outputs[1])
    pool  = CachingPool(workers())
    if progressbar(outputs[1])
        lck = ReentrantLock()
        prg = Progress(length(queue))
        increment() = (lock(lck); next!(prg); unlock(lck))
    end
    @sync for worker in workers()
        @async while length(queue) > 0
            indices = shift!(queue)
            remotecall_fetch(just_do_it!, pool, f, outputs, args, indices)
            progressbar(outputs[1]) && increment()
        end
    end
    output
end

function just_do_it!(f, output, args, indices) # (c) Nike
    output[indices...] = f(getindex.(args, indices...)...)
end

function just_do_it!(f, outputs::Tuple, args, indices) # (c) Nike
    result = f(getindex.(args, indices...)...)
    for idx = 1:length(result)
        output[idx][indices...] = result[idx]
    end
end

