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

function Base.broadcast!(f, output::BlockMatrix, args...)
    flatten(x) = reshape(x, length(x))
    queue = flatten(collect(Iterators.product(indices(output)...)))
    pool  = CachingPool(workers())
    lck = ReentrantLock()
    prg = Progress(length(queue))
    increment() = (lock(lck); next!(prg); unlock(lck))
    @sync for worker in workers()
        @async while length(queue) > 0
            indices = shift!(queue)
            remotecall_fetch(do_on_remote!, pool, f, output, args, indices)
            increment()
        end
    end
    output
end

function do_on_remote!(f, output, args, indices)
    output[indices...] = f(getindex.(args, indices...)...)
end

