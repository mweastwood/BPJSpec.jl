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

struct Workers
    dict :: Dict{String, Vector{Int}}
end

function categorize_workers()
    futures = [remotecall(() -> chomp(readstring(`hostname`)), worker) for worker in workers()]
    hierarchy = Dict{String, Vector{Int}}()
    for (future, worker) in zip(futures, workers())
        hostname = fetch(future)
        if haskey(hierarchy, hostname)
            push!(hierarchy[hostname], worker)
        else
            hierarchy[hostname] = [worker]
        end
    end
    Workers(hierarchy)
end

function Base.show(io::IO, workers::Workers)
    hosts = collect(keys(workers.dict))
    sort!(hosts)
    println("| Workers")
    println("|---------")
    for host in hosts
        ids = workers.dict[host]
        @printf(io, "| %7s : %s\n", host, join(string.(ids), ", "))
    end
end

#"""
#    @distribute
#
#Use this macro in combination with `@remote` to distribute the work
#from a for-loop across all available workers.
#
## Example
#
#    @distribute for idx = 1:100
#        result = @remote do_something_on_a_worker()
#        do_something_locally()
#    end
#"""
#macro distribute(for_loop)
#    idx = for_loop.args[1].args[1]
#    lower_limit = for_loop.args[1].args[2].args[1]
#    upper_limit = for_loop.args[1].args[2].args[2]
#    inner_block = for_loop.args[2]
#    expr = quote
#        $idx = $lower_limit
#        nextidx() = (myidx = $idx; $idx += 1; myidx)
#        p = Progress($upper_limit - $lower_limit + 1, "Progress: ")
#        l = ReentrantLock()
#        increment_progress() = (lock(l); next!(p); unlock(l))
#        @sync for worker in workers()
#            @async while true
#                $idx = nextidx()
#                $idx ≤ $upper_limit || break
#                $inner_block
#                increment_progress()
#            end
#        end
#    end
#    esc(expr)
#end
#
#macro remote(function_call)
#    func = function_call.args[1]
#    args = function_call.args[2:end]
#    expr = quote
#        remotecall_fetch(worker, $func, $(args...))
#    end
#    esc(expr)
#end

