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

"Wraps a type `T` such that a progress bar is displayed when broadcasting into `T`."
struct ProgressBar{T}
    data :: T
end

unwrap(x::ProgressBar) = x.data
unwrap(x) = x

const BroadcastTypes = Union{ProgressBar, AbstractBlockMatrix}

function Base.broadcast!(f, output::BroadcastTypes, args...)
    internal_broadcast!(f, output, unwrap.(args)...)
end

function internal_broadcast!(f, output::ProgressBar, args...)
    internal_broadcast!(f, unwrap(output), args..., progress=true)
end

function internal_broadcast!(f, output, args...; progress=false)
    distribute = distribute_write(output) && all(distribute_read.(args))
    if distribute
        return distributed_broadcast!(f, output, args, progress=progress)
    else
        return local_broadcast!(f, output, args, progress=progress)
    end
end

function local_broadcast!(f, output, args; progress=false)
    queue = collect(indices(output))
    progress && (prg = Progress(length(queue)))
    for indices in queue
        just_do_it!(f, output, args, indices)
        progress && next!(prg)
    end
    output
end

function distributed_broadcast!(f, output, args; progress=false)
    queue = collect(indices(output))
    pool  = CachingPool(workers())
    if progress
        lck = ReentrantLock()
        prg = Progress(length(queue))
        increment() = (lock(lck); next!(prg); unlock(lck))
    end
    @sync for worker in workers()
        @async while length(queue) > 0
            indices = shift!(queue)
            remotecall_fetch(just_do_it!, pool, f, output, args, indices)
            progress && increment()
        end
    end
    output
end

function multi_broadcast!(f, outputs, args; progress=false)
    queue = collect(indices(outputs[1]))
    pool  = CachingPool(workers())
    if progress
        lck = ReentrantLock()
        prg = Progress(length(queue))
        increment() = (lock(lck); next!(prg); unlock(lck))
    end
    @sync for worker in workers()
        @async while length(queue) > 0
            indices = shift!(queue)
            remotecall_fetch(just_do_it!, pool, f, outputs, args, indices)
            progress && increment()
        end
    end
    outputs
end

function just_do_it!(f, output, args, indices) # (c) Nike
    output[indices...] = f(getindex.(args, indices...)...)
end

function just_do_it!(f, outputs::Tuple, args, indices) # (c) Nike
    local results
    try
        results = f(getindex.(args, indices...)...)
    catch exception
        @show indices
        rethrow(exception)
    end
    for idx = 1:length(results)
        outputs[idx][indices...] = results[idx]
    end
end

