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

abstract type TransferMatrix end

struct HierarchicalTransferMatrix <: TransferMatrix
    path :: String
    metadata :: Metadata
    beam :: Function
    function HierarchicalTransferMatrix(path, metadata, beam)
        isdir(path) || mkdir(path)
        save(joinpath(path, "METADATA.jld2"), "metadata", metadata, "beam", beam)
        new(path, metadata, beam)
    end
end

function HierarchicalTransferMatrix(path)
    metadata, beam = load(joinpath(path, "METADATA.jld2"), "metadata", "beam")
    HierarchicalTransferMatrix(path, metadata, beam)
end

function compute!(transfermatrix::HierarchicalTransferMatrix)
    println("")
    println("| Starting transfer matrix calculation *")
    println("|---------")
    println("| ($(now()))")
    println("")

    workers = categorize_workers()
    println(workers)

    hierarchy = compute_baseline_hierarchy(transfermatrix.metadata)
    println(hierarchy)
    save(joinpath(transfermatrix.path, "HIERARCHY.jld2"), "hierarchy", hierarchy)

    #for ν in transfermatrix.metadata.frequencies
        ν = transfermatrix.metadata.frequencies[55]
        compute_one_frequency!(transfermatrix, workers, hierarchy, ν)
    #end
end

function compute_one_frequency!(transfermatrix::HierarchicalTransferMatrix,
                                workers, hierarchy, ν)
    metadata = transfermatrix.metadata

    for idx = 1:length(hierarchy.divisions)-1
        lmax = hierarchy.divisions[idx+1]
        baselines = transfermatrix.metadata.baselines[hierarchy.baselines[idx]]
        compute_baseline_group_one_frequency!(transfermatrix, workers.dict["astm11"],
                                              baselines, lmax, ν)
    end
end

function compute_baseline_group_one_frequency!(transfermatrix::HierarchicalTransferMatrix,
                                               workers, baselines, lmax, ν)
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
    beam = create_beam_map(transfermatrix.beam, transfermatrix.metadata, (lmax+1, 2mmax+1))
    rhat = unit_vectors(beam)
    plan = plan_sht(lmax, mmax, size(rhat))

    queue  = collect(1:length(baselines))
    blocks = [zeros(Complex128, two(m)*length(baselines), lmax-m+1) for m = 0:mmax]

    function just_do_it(α)
        real_coeff, imag_coeff = fringe_pattern(baselines[α], phase_center, beam, rhat, plan, ν)
    end

    lck = ReentrantLock()
    prg = Progress(length(queue))
    increment() = (lock(lck); next!(prg); unlock(lck))

    @sync for worker in workers
        @async while length(queue) > 0
            α = pop!(queue)
            real_coeff, imag_coeff = remotecall_fetch(just_do_it, pool, α)
            write_to_blocks!(blocks, real_coeff, imag_coeff, lmax, mmax, α)
            increment()
        end
    end

    path = joinpath(transfermatrix.path, @sprintf("%.3fMHz", ustrip(uconvert(u"MHz", ν))))
    isdir(path) || mkdir(path)
    path = joinpath(path, @sprintf("lmax=%04d", lmax))
    isdir(path) || mkdir(path)

    prg = Progress(mmax+1)
    for m = 0:mmax
        save(joinpath(path, @sprintf("m=%04d.jld2", m)), "block", blocks[m+1], compress=true)
        next!(prg)
    end

    blocks
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
function fringe_pattern(baseline, phase_center, beam, rhat, plan, ν)
    λ = u"c" / ν
    real_fringe, imag_fringe = plane_wave(rhat, baseline, phase_center, λ)
    real_coeff = plan * Map(real_fringe .* beam)
    imag_coeff = plan * Map(imag_fringe .* beam)
    real_coeff, imag_coeff
end

function plane_wave(rhat, baseline, phase_center, λ)
    real_part = similar(rhat, Float64)
    imag_part = similar(rhat, Float64)
    two_π = 2π
    for idx in eachindex(rhat)
        ϕ = uconvert(u"rad", two_π*dot(rhat[idx] - phase_center, baseline)/λ)
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

"Creat an image of the beam model."
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





#function generate_transfermatrix_onechannel!(transfermatrix, meta, beam, variables, ν)
#    lmax = transfermatrix.lmax
#    mmax = transfermatrix.mmax
#    # Memory map all the blocks on the master process to avoid having to
#    # open/close the files multiple times and to avoid having to read the
#    # entire matrix at once.
#    info("Running new version!")
#    info("Memory mapping files")
#    blocks = IOStream[]
#    #blocks = Matrix{Complex128}[]
#    for m = 0:mmax
#        directory = directory_name(m, ν, mmax+1)
#        directory = joinpath(transfermatrix.path, directory)
#        isdir(directory) || mkdir(directory)
#        filename = block_filename(m, ν)
#        block = open(joinpath(directory, filename), "w")
#        # Write the size of the matrix block to the start of
#        # the file because in general we don't know how many
#        # rows will be in each block of the matrix.
#        #
#        # Also note that we are storing the transpose of each
#        # block in order to make all the disk writes sequential.
#        sz = (lmax-m+1, two(m)*Nbase(meta))
#        write(block, sz[1], sz[2])
#        push!(blocks, block)
#        #open(joinpath(directory, filename), "w+") do file
#        #    # note that we store the transpose of the transfer matrix blocks to make
#        #    # all the disk writes sequential
#        #    sz = (lmax-m+1, two(m)*Nbase(meta))
#        #    write(file, sz[1], sz[2])
#        #    block = Mmap.mmap(file, Matrix{Complex128}, sz)
#        #    push!(blocks, block)
#        #end
#        #open(joinpath(directory, filename), "r+") do file
#        #    sz1 = read(file, Int)
#        #    sz2 = read(file, Int)
#        #    sz = (sz1, sz2)
#        #    block = Mmap.mmap(file, Matrix{Complex128}, sz)
#        #    push!(blocks, block)
#        #end
#    end
#    info("Beginning the computation")
#    idx = 1
#    #idx = 1500
#    nextidx() = (myidx = idx; idx += 1; myidx)
#    p = Progress(Nbase(meta) - idx + 1, "Progress: ")
#    l = ReentrantLock()
#    increment_progress() = (lock(l); next!(p); unlock(l))
#    @sync for worker in workers()
#        @async begin
#            input = RemoteChannel()
#            output_realfringe = RemoteChannel()
#            output_imagfringe = RemoteChannel()
#            remotecall(transfermatrix_worker_loop, worker,
#                       input, output_realfringe, output_imagfringe, beam, variables, ν)
#            while true
#                α = nextidx()
#                α ≤ Nbase(meta) || break
#                put!(input, α)
#                realfringe = take!(output_realfringe)
#                imagfringe = take!(output_imagfringe)
#                pack!(blocks, realfringe, imagfringe, lmax, mmax, α)
#                increment_progress()
#            end
#        end
#    end
#    for block in blocks
#        close(block)
#    end
#end
#
#function transfermatrix_worker_loop(input, output_realfringe, output_imagfringe, beam, variables, ν)
#    while true
#        α = take!(input)
#        realfringe, imagfringe = fringes(beam, variables, ν, α)
#        put!(output_realfringe, realfringe)
#        put!(output_imagfringe, imagfringe)
#    end
#end
#
#"""
#    fringes(beam, variables, ν, α)
#
#Generate the spherical harmonic expansion of the fringe pattern on the sky.
#
#Note that because the Healpix library assumes you are asking for the coefficients
#of a real field, there must be one set of coefficients for the real part of
#the fringe pattern and one set of coefficients for the imaginary part of the
#fringe pattern.
#"""
#function fringes(beam, variables, ν, α)
#    λ = c / ν
#    u = variables.u[α] / λ
#    v = variables.v[α] / λ
#    w = variables.w[α] / λ
#    realmap, imagmap = planewave(u, v, w, variables.x, variables.y, variables.z, variables.phase_center)
#    realfringe = map2alm(beam .* realmap, variables.lmax, variables.mmax, iterations=2)
#    imagfringe = map2alm(beam .* imagmap, variables.lmax, variables.mmax, iterations=2)
#    realfringe, imagfringe
#end
#
#"""
#    pack!(blocks, realfringe, imagfringe, lmax, mmax, α)
#
#Having calculated the spherical harmonic expansion of the fringe pattern,
#pack those numbers into the transfer matrix.
#"""
#function pack!(blocks, realfringe, imagfringe, lmax, mmax, α)
#    # Note that all the conjugations in this function come about because
#    # Shaw et al. 2014, 2015 expand the fringe pattern in terms of the
#    # spherical harmonic conjugates while we've expanded the fringe pattern
#    # in terms of the spherical harmonics.
#    #for l = 0:lmax
#    #    blocks[1][l+1,α] = conj(realfringe[l,0]) + 1im*conj(imagfringe[l,0])
#    #end
#    #for m = 1:mmax
#    #    block = blocks[m+1]
#    #    α1 = 2α-1 # positive m
#    #    for l = m:lmax
#    #        block[l-m+1,α1] = conj(realfringe[l,m]) + 1im*conj(imagfringe[l,m])
#    #    end
#    #    α2 = 2α-0 # negative m
#    #    for l = m:lmax
#    #        block[l-m+1,α2] = conj(realfringe[l,m]) - 1im*conj(imagfringe[l,m])
#    #    end
#    #end
#    offset = 2sizeof(Int) + (α-1)*(lmax+1)*sizeof(Complex128)
#    output = Complex128[conj(realfringe[l,0]) + 1im*conj(imagfringe[l,0]) for l = 0:lmax]
#    seek(blocks[1], offset)
#    write(blocks[1], output)
#    for m = 1:mmax
#        offset = 2sizeof(Int) + 2*(α-1)*(lmax-m+1)*sizeof(Complex128)
#        output1 = Complex128[conj(realfringe[l,m]) + 1im*conj(imagfringe[l,m]) for l = m:lmax] # positive m
#        output2 = Complex128[conj(realfringe[l,m]) - 1im*conj(imagfringe[l,m]) for l = m:lmax] # negative m
#        seek(blocks[m+1], offset)
#        write(blocks[m+1], output1, output2)
#    end
#end
#
#function setindex!(transfermatrix::TransferMatrix, block, m, channel)
#    ν = transfermatrix.frequencies[channel]
#    directory = directory_name(m, ν, transfermatrix.mmax+1)
#    filename = block_filename(m, ν)
#    open(joinpath(transfermatrix.path, directory, filename), "w") do file
#        write(file, size(block, 2), size(block, 1), block.')
#    end
#    block
#end
#
#function getindex(transfermatrix::TransferMatrix, m, channel)
#    local block
#    ν = transfermatrix.frequencies[channel]
#    directory = directory_name(m, ν, transfermatrix.mmax+1)
#    filename = block_filename(m, ν)
#    open(joinpath(transfermatrix.path, directory, filename), "r") do file
#        sz = tuple(read(file, Int, 2)...)
#        block = read(file, Complex128, sz)
#    end
#    block.'
#end
#
#
#doc"""
#    preserve_singular_values(B::TransferMatrix)
#
#Construct a matrix that projects the $m$-modes onto a lower dimensional
#space while preserving all the singular values of the transfer matrix.
#
#Multiplying by this matrix will compress the data, make the transfer
#matrix square, and leave the information about the sky untouched.
#"""
#function preserve_singular_values(B::TransferMatrix)
#    N = length(B.blocks)
#    blocks = Array{MatrixBlock}(N)
#    for i = 1:N
#        U,σ,V = svd(B.blocks[i])
#        blocks[i] = MatrixBlock(U')
#    end
#    Blocks(blocks)
#end
#

