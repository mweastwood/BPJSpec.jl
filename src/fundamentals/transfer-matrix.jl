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

doc"""
    struct TransferMatrix

This singleton type represents the transfer matrix of an interferometer. This matrix effectively
describes how an interferometer responds to the sky, including the antenna primary beam, bandpass,
and baseline distribution.

This matrix is hierarchical in the sense that we save on some computational and storage requirements
by separating long baselines from short baselines.

**Usage:**

**See also:** [`MModes`](@ref), [`NoiseCovarianceMatrix`](@ref), [`MFBlockMatrix`](@ref)
"""
struct TransferMatrix end

function create(::Type{TransferMatrix}, path::String, metadata::Metadata, beam;
                lmax=-1, compute=true, rm=false, progress=false)
    hierarchy = Hierarchy(metadata, lmax=lmax)
    storage   = HierarchicalStorage(path, hierarchy)
    output = create(MFBlockMatrix, storage, hierarchy.lmax,
                    metadata.frequencies, metadata.bandwidth, rm=rm)
    compute && compute!(TransferMatrix, output, metadata, beam, progress=progress)
    output
end

function compute!(::Type{TransferMatrix}, matrix::MFBlockMatrix{HierarchicalStorage},
                  metadata::Metadata, beam; progress=false)
    if progress
        println("")
        println("| Starting transfer matrix calculation")
        println("|---------")
        println("| ($(now()))")
        println("")
    end

    workers = categorize_workers()

    if progress
        println(workers)
        println(matrix.storage.hierarchy)
    end

    queue = collect(1:length(metadata.frequencies))
    if progress
        lck = ReentrantLock()
        prg = Progress(length(queue))
        increment() = (lock(lck); next!(prg); unlock(lck))
    end

    @sync for worker in leaders(workers)
        @async while length(queue) > 0
            β = shift!(queue)
            remotecall_fetch(compute_one_frequency!, worker, matrix, workers, metadata, beam, β)
            progress && increment()
        end
    end
end

function compute_one_frequency!(matrix, workers, metadata, beam, β)
    hierarchy = matrix.storage.hierarchy

    my_machine   = chomp(readstring(`hostname`))
    subordinates = copy(workers.dict[my_machine])
    if length(subordinates) > 1
        # make sure this process isn't in the worker pool
        deleteat!(subordinates, subordinates .== myid())
    end

    for idx = 1:length(hierarchy.divisions)-1
        lmax = hierarchy.divisions[idx+1]
        baselines = metadata.baselines[hierarchy.baselines[idx]]
        blocks = compute_baseline_group_one_frequency!(matrix, subordinates,
                                                       metadata, beam, baselines, lmax, β)
        Base.resize!(blocks, 0)
        finalize(blocks)
        gc(); gc() # please please please garbage collect `blocks`
    end
end

function compute_baseline_group_one_frequency!(matrix, subordinates, metadata,
                                               beam, baselines, lmax, β)
    pool = CachingPool(subordinates)

    # "... but in this world nothing can be said to be certain, except death
    #  and taxes and lmax=mmax"
    #   - Benjamin Franklin, 1789
    mmax = lmax
    ν = metadata.frequencies[β]
    phase_center = metadata.phase_center
    beam_map = create_beam_map(beam, metadata, (lmax+1, 2mmax+1))
    rhat = unit_vectors(beam_map)
    plan = FastTransformsWrapper.plan_sht(lmax, mmax, size(rhat))

    queue  = collect(1:length(baselines))
    blocks = [zeros(Complex128, two(m)*length(baselines), lmax-m+1) for m = 0:mmax]

    function just_do_it(α)
        real_coeff, imag_coeff = fringe_pattern(baselines[α], phase_center, beam_map, rhat, plan, ν)
    end

    @sync for subordinate in subordinates
        @async while length(queue) > 0
            α = pop!(queue)
            real_coeff, imag_coeff = remotecall_fetch(just_do_it, pool, α)
            fix_scaling!(real_coeff, imag_coeff, ν)
            write_to_blocks!(blocks, real_coeff, imag_coeff, lmax, mmax, α)
        end
    end

    for m = 0:mmax
        matrix.storage[lmax, m, β] = blocks[m+1]
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
    real_coeff = plan * FastTransformsWrapper.Map(real_fringe .* beam_map)
    imag_coeff = plan * FastTransformsWrapper.Map(imag_fringe .* beam_map)
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
    FastTransformsWrapper.Map(real_part), FastTransformsWrapper.Map(imag_part)
end

"Compute the unit vector to each point on the sky."
function unit_vectors(map)
    rhat = Matrix{Direction}(size(map))
    for jdx = 1:size(map, 2), idx = 1:size(map, 1)
        rhat[idx, jdx] = FastTransformsWrapper.index2vector(map, idx, jdx)
    end
    rhat
end

"Create an image of the beam model."
function create_beam_map(f, metadata, size)
    zenith = Direction(metadata.position)
    north  = gram_schmidt(Direction(dir"ITRF", 0, 0, 1), zenith)
    east   = cross(north, zenith)

    map = FastTransformsWrapper.Map(zeros(size))
    for jdx = 1:size[2], idx = 1:size[1]
        vec = FastTransformsWrapper.index2vector(map, idx, jdx)
        x = dot(vec, east)
        y = dot(vec, north)
        z = dot(vec, zenith)
        elevation = asin(clamp(z, -1, 1))
        azimuth   = atan2(x, y)
        map[idx, jdx] = f(azimuth, elevation)
    end
    map
end

