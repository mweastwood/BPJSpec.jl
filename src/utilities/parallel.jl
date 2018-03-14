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
    println(io, "| Workers")
    println(io, "|---------")
    for host in hosts
        ids = workers.dict[host]
        @printf(io, "| %7s : %s\n", host, join(string.(ids), ", "))
    end
end

"Return a list of one worker per machine."
leaders(workers::Workers) = first.(collect(values(workers.dict)))

