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

abstract type TransferMatrix <: BlockMatrix end

struct HierarchicalTransferMatrix <: TransferMatrix
    path     :: String
    metadata :: Metadata
    function HierarchicalTransferMatrix(path, metadata)
        isdir(path) || mkpath(path)
        save(joinpath(path, "METADATA.jld2"), "metadata", metadata)
        new(path, metadata)
    end
end

function HierarchicalTransferMatrix(path)
    metadata = load(joinpath(path, "METADATA.jld2"), "metadata")
    HierarchicalTransferMatrix(path, metadata)
end

function compute!(transfermatrix::HierarchicalTransferMatrix, beam;
                  lmax=maximum(maximum_multipole_moment(transfermatrix.metadata)))
    println("")
    println("| Starting transfer matrix calculation")
    println("|---------")
    println("| ($(now()))")
    println("")

    workers = categorize_workers()
    println(workers)

    hierarchy = compute_baseline_hierarchy(transfermatrix.metadata, lmax)
    println(hierarchy)
    save(joinpath(transfermatrix.path, "HIERARCHY.jld2"), "hierarchy", hierarchy)

    for ν in transfermatrix.metadata.frequencies[11:end]
        @show ν
        @time compute_one_frequency!(transfermatrix, workers, hierarchy, beam, ν)
    end
end

function compute_one_frequency!(transfermatrix::HierarchicalTransferMatrix,
                                workers, hierarchy, beam, ν)
    metadata = transfermatrix.metadata

    for idx = 1:length(hierarchy.divisions)-1
        @time begin
            lmax = hierarchy.divisions[idx+1]
            baselines = transfermatrix.metadata.baselines[hierarchy.baselines[idx]]
            blocks = compute_baseline_group_one_frequency!(transfermatrix, workers.dict["astm11"],
                                                           beam, baselines, lmax, ν)
            resize!(blocks, 0)
            finalize(blocks)
            gc(); gc() # please please please garbage collect `blocks`
        end
    end
end

function compute_baseline_group_one_frequency!(transfermatrix::HierarchicalTransferMatrix,
                                               workers, beam, baselines, lmax, ν)
    if length(workers) > 1
        # make sure this process isn't in the worker pool
        deleteat!(workers, workers .== myid())
    end
    pool = CachingPool(workers)

    # "... but in this world nothing can be said to be certain, except death
    #  and taxes and lmax=mmax"
    #   - Benjamin Franklin, 1789
    mmax = lmax
    phase_center = transfermatrix.metadata.phase_center
    beam_map = create_beam_map(beam, transfermatrix.metadata, (lmax+1, 2mmax+1))
    rhat = unit_vectors(beam_map)
    plan = plan_sht(lmax, mmax, size(rhat))

    queue  = collect(1:length(baselines))
    blocks = [zeros(Complex128, two(m)*length(baselines), lmax-m+1) for m = 0:mmax]

    function just_do_it(α)
        real_coeff, imag_coeff = fringe_pattern(baselines[α], phase_center, beam_map, rhat, plan, ν)
    end

    lck = ReentrantLock()
    prg = Progress(length(queue))
    increment() = (lock(lck); next!(prg); unlock(lck))

    @sync for worker in workers
        @async while length(queue) > 0
            α = pop!(queue)
            real_coeff, imag_coeff = remotecall_fetch(just_do_it, pool, α)
            fix_scaling!(real_coeff, imag_coeff, ν)
            write_to_blocks!(blocks, real_coeff, imag_coeff, lmax, mmax, α)
            increment()
        end
    end

    path = joinpath(transfermatrix.path, @sprintf("%.3fMHz", ustrip(uconvert(u"MHz", ν))))
    isdir(path) || mkdir(path)

    prg = Progress(mmax+1)
    jldopen(joinpath(path, @sprintf("lmax=%04d.jld2", lmax)), "w") do file
        for m = 0:mmax
            file[@sprintf("%04d", m)] = blocks[m+1]
            next!(prg)
        end
    end

    blocks
end

function fix_scaling!(real_coeff, imag_coeff, ν)
    # Our m-modes are in units of Jy, but our alm are in units of K. Here we apply the scaling
    # factor to the transfer matrix that makes this work with the right units.

    # This is the conversion factor I have been using to convert my alm into units of K. We'll apply
    # the inverse here to the transfer matrix so that this conversion factor is no longer necessary.
    factor = ustrip(uconvert(u"K", u"Jy * c^2/(2*k)"/ν^2))
    real_coeff.matrix ./= factor
    imag_coeff.matrix ./= factor
end

function write_to_blocks!(blocks, real_coeff, imag_coeff, lmax, mmax, α)
    # m = 0
    block = blocks[1]
    for l = 0:lmax
        block[α, l+1] = conj(real_coeff[l, 0]) + 1im*conj(imag_coeff[l, 0])
    end
    # m > 0
    for m = 1:mmax
        block = blocks[m+1]
        α1 = 2α-1 # positive m
        α2 = 2α-0 # negative m
        for l = m:lmax
            block[α1, l-m+1] = conj(real_coeff[l, m]) + 1im*conj(imag_coeff[l, m])
            block[α2, l-m+1] = conj(real_coeff[l, m]) - 1im*conj(imag_coeff[l, m])
        end
    end
end

"Compute the spherical harmonic transform of the fringe pattern for the given baseline."
function fringe_pattern(baseline, phase_center, beam_map, rhat, plan, ν)
    λ = u"c" / ν
    real_fringe, imag_fringe = plane_wave(rhat, baseline, phase_center, λ)
    real_coeff = plan * Map(real_fringe .* beam_map)
    imag_coeff = plan * Map(imag_fringe .* beam_map)
    real_coeff, imag_coeff
end

function plane_wave(rhat, baseline, phase_center, λ)
    real_part = similar(rhat, Float64)
    imag_part = similar(rhat, Float64)
    two_π = 2π
    δϕ = two_π*dot(phase_center, baseline)/λ
    for idx in eachindex(rhat)
        ϕ = uconvert(u"rad", two_π*dot(rhat[idx], baseline)/λ - δϕ)
        real_part[idx] = cos(ϕ)
        imag_part[idx] = sin(ϕ)
    end
    Map(real_part), Map(imag_part)
end

"Compute the unit vector to each point on the sky."
function unit_vectors(map)
    rhat = Matrix{Direction}(size(map))
    for jdx = 1:size(map, 2), idx = 1:size(map, 1)
        rhat[idx, jdx] = index2vector(map, idx, jdx)
    end
    rhat
end

"Create an image of the beam model."
function create_beam_map(f, metadata, size)
    zenith = Direction(metadata.position)
    north  = gram_schmidt(Direction(dir"ITRF", 0, 0, 1), zenith)
    east   = cross(north, zenith)

    map = BPJSpec.Map(zeros(size))
    for jdx = 1:size[2], idx = 1:size[1]
        vec = index2vector(map, idx, jdx)
        x = dot(vec, east)
        y = dot(vec, north)
        z = dot(vec, zenith)
        elevation = asin(clamp(z, -1, 1))
        azimuth   = atan2(x, y)
        map[idx, jdx] = f(azimuth, elevation)
    end
    map
end

function Base.getindex(transfermatrix::HierarchicalTransferMatrix, m, β)
    ν = transfermatrix.metadata.frequencies[β]

    # load each hierarchical component of the transfer matrix
    hierarchy = load(joinpath(transfermatrix.path, "HIERARCHY.jld2"), "hierarchy")
    path = joinpath(transfermatrix.path, @sprintf("%.3fMHz", ustrip(uconvert(u"MHz", ν))))
    blocks = Matrix{Complex128}[]
    for idx = 1:length(hierarchy.divisions)-1
        lmax = hierarchy.divisions[idx+1]
        m > lmax && continue
        filename   = @sprintf("lmax=%04d.jld2", lmax)
        objectname = @sprintf("%04d", m)
        push!(blocks, load(joinpath(path, filename), objectname))
    end

    # stitch the components together into a single matrix
    output = zeros(Complex128, sum(size(block, 1) for block in blocks),
                   maximum(size(block, 2) for block in blocks))
    offset = 1
    for block in blocks
        range1 = offset:offset+size(block, 1)-1
        range2 = 1:size(block, 2)
        output_view = @view output[range1, range2]
        copy!(output_view, block)
        offset += size(block, 1)
    end
    output
end

function getlmax(transfermatrix::HierarchicalTransferMatrix)
    # TODO: I'd like to have this stored as a field...
    hierarchy = load(joinpath(transfermatrix.path, "HIERARCHY.jld2"), "hierarchy")
    hierarchy.divisions[end]
end

"Get the baseline permutation vector for the given value of m."
function baseline_permutation(transfermatrix::HierarchicalTransferMatrix, m)
    hierarchy = load(joinpath(transfermatrix.path, "HIERARCHY.jld2"), "hierarchy")
    indices = Int[]
    for idx = 1:length(hierarchy.divisions)-1
        lmax = hierarchy.divisions[idx+1]
        m > lmax && continue
        if m == 0
            append!(indices, hierarchy.baselines[idx])
        else
            for baseline in hierarchy.baselines[idx]
                push!(indices, 2*baseline-1) # positive m
                push!(indices, 2*baseline-0) # negative m
            end
        end
    end
    indices
end

