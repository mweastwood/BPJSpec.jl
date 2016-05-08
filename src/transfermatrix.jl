# Copyright (c) 2015 Michael Eastwood
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

function TransferMatrix(path)
    local lmax, mmax, frequencies
    open(joinpath(path, "METADATA"), "r") do file
        lmax = read(file, Int)
        mmax = read(file, Int)
        len  = read(file, Int)
        frequencies = read(file, Float64, len)
    end
    TransferMatrix(path, lmax, mmax, frequencies)
end

Nfreq(matrix::TransferMatrix) = length(matrix.frequencies)

"""
    TransferMatrix(path, meta::Metadata, lmax, mmax, nside)

Generate the transfer matrix. This will take some time.
"""
function TransferMatrix(path, meta::Metadata, lmax, mmax, nside)
    transfermatrix = TransferMatrix(path, lmax, mmax, meta.channels)
    variables = TransferMatrixVariables(meta, lmax, mmax, nside)
    generate_transfermatrix!(transfermatrix, meta, variables)
    transfermatrix
end

function TransferMatrix(path, meta::Metadata, beam::HealpixMap, lmax, mmax, nside)
    transfermatrix = TransferMatrix(path, lmax, mmax, meta.channels)
    variables = TransferMatrixVariables(meta, lmax, mmax, nside)
    generate_transfermatrix!(transfermatrix, meta, beam, variables)
    transfermatrix
end

"""
    TransferMatrixVariables

This type just holds some variables we would like to precompute
and use while generating the transfer matrix.

# Fields

* `lmax` and `mmax` specify the space of spherical harmonics to use
* `x`, `y`, and `z` give the unit vector to each pixel in the Healpix map
* `u`, `v`, and `w` give each baseline vector in ITRF coordinates (meters)
* `phase_center` the ITRF direction to the phase center
"""
immutable TransferMatrixVariables
    lmax :: Int
    mmax :: Int
    x :: HealpixMap
    y :: HealpixMap
    z :: HealpixMap
    u :: Vector{Float64}
    v :: Vector{Float64}
    w :: Vector{Float64}
    phase_center :: Direction
end

function TransferMatrixVariables(meta, lmax, mmax, nside)
    x = HealpixMap(Float64, nside)
    y = HealpixMap(Float64, nside)
    z = HealpixMap(Float64, nside)
    for pix = 1:length(x)
        vec = LibHealpix.pix2vec_ring(nside, pix)
        x[pix] = vec[1]
        y[pix] = vec[2]
        z[pix] = vec[3]
    end

    nbase = Nbase(meta)
    u = zeros(nbase)
    v = zeros(nbase)
    w = zeros(nbase)
    for α = 1:nbase
        antenna1 = meta.antennas[meta.baselines[α].antenna1]
        antenna2 = meta.antennas[meta.baselines[α].antenna2]
        u[α] = antenna1.position.x - antenna2.position.x
        v[α] = antenna1.position.y - antenna2.position.y
        w[α] = antenna1.position.z - antenna2.position.z
    end

    frame = TTCal.reference_frame(meta)
    phase_center = measure(frame, meta.phase_center, dir"ITRF")

    TransferMatrixVariables(lmax, mmax, x, y, z, u, v, w, phase_center)
end

function generate_transfermatrix!(transfermatrix, meta, variables)
    for ν in transfermatrix.frequencies
        beam = beam_map(meta, ν)
        generate_transfermatrix_onechannel!(transfermatrix, meta, beam, variables, ν)
    end
end

function generate_transfermatrix!(transfermatrix, meta, beam, variables)
    for ν in transfermatrix.frequencies
        generate_transfermatrix_onechannel!(transfermatrix, meta, beam, variables, ν)
    end
end

function generate_transfermatrix_onechannel!(transfermatrix, meta, beam, variables, ν)
    lmax = transfermatrix.lmax
    mmax = transfermatrix.mmax
    # Memory map all the blocks on the master process to avoid having to
    # open/close the files multiple times and to avoid having to read the
    # entire matrix at once.
    info("Memory mapping files")
    blocks = Matrix{Complex128}[]
    for m = 0:transfermatrix.mmax
        filename = block_filename(m, ν)
        open(joinpath(transfermatrix.path, filename), "w+") do file
            # note that we store the transpose of the transfer matrix blocks to make
            # all the disk writes sequential
            sz = (lmax-m+1, two(m)*Nbase(meta))
            write(file, sz[1], sz[2])
            block = Mmap.mmap(file, Matrix{Complex128}, sz)
            push!(blocks, block)
        end
    end
    info("Beginning the computation")
    @distribute for α = 1:Nbase(meta)
        realfringe, imagfringe = @remote fringes(beam, variables, ν, α)
        pack!(blocks, realfringe, imagfringe, lmax, mmax, α)
    end
end

"""
    beam_map(meta::Metadata, frequency)

Create a Healpix map of the Stokes I beam.
"""
function beam_map(meta::Metadata, frequency)
    frame  = TTCal.reference_frame(meta)
    position = measure(frame, TTCal.position(meta), pos"ITRF")
    zenith = normalize!([position.x, position.y, position.z])
    north  = gramschmidt([0.0, 0.0, 1.0], zenith)
    east   = cross(north, zenith)
    map = HealpixMap(zeros(nside2npix(512)))
    for i = 1:length(map)
        vec = LibHealpix.pix2vec_ring(512, i)
        el = π/2 - angle_between(vec, zenith)
        x  = dot(vec, east)
        y  = dot(vec, north)
        az = atan2(x, y)
        if el < 0
            map[i] = 0
        else
            J = meta.beam(frequency, az, el)
            M = MuellerMatrix(J)
            map[i] = M.mat[1,1]
        end
    end
    map
end

"""
    planewave(u, v, w, x, y, z, phase_center)

Compute the fringe pattern over a Healpix image.

```math
exp(2 \pi i (ux+vy+wz)
```
"""
function planewave(u, v, w, x, y, z, phase_center)
    realmap = HealpixMap(Float64, nside(x))
    imagmap = HealpixMap(Float64, nside(x))
    for idx = 1:length(realmap)
        δx = x[idx] - phase_center.x
        δy = y[idx] - phase_center.y
        δz = z[idx] - phase_center.z
        ϕ = 2π*(u*δx + v*δy + w*δz)
        realmap[idx] = cos(ϕ)
        imagmap[idx] = sin(ϕ)
    end
    realmap, imagmap
end

"""
    fringes(beam, variables, ν, α)

Generate the spherical harmonic expansion of the fringe pattern on the sky.

Note that because the Healpix library assumes you are asking for the coefficients
of a real field, there must be one set of coefficients for the real part of
the fringe pattern and one set of coefficients for the imaginary part of the
fringe pattern.
"""
function fringes(beam, variables, ν, α)
    λ = c / ν
    u = variables.u[α] / λ
    v = variables.v[α] / λ
    w = variables.w[α] / λ
    realmap, imagmap = planewave(u, v, w, variables.x, variables.y, variables.z, variables.phase_center)
    realfringe = map2alm(beam .* realmap, variables.lmax, variables.mmax, iterations=5)
    imagfringe = map2alm(beam .* imagmap, variables.lmax, variables.mmax, iterations=5)
    realfringe, imagfringe
end

"""
    pack!(blocks, realfringe, imagfringe, lmax, mmax, α)

Having calculated the spherical harmonic expansion of the fringe pattern,
pack those numbers into the transfer matrix.
"""
function pack!(blocks, realfringe, imagfringe, lmax, mmax, α)
    # Note that all the conjugations in this function come about because
    # Shaw et al. 2014, 2015 expand the fringe pattern in terms of the
    # spherical harmonic conjugates while we've expanded the fringe pattern
    # in terms of the spherical harmonics.
    for l = 0:lmax
        blocks[1][l+1,α] = conj(realfringe[l,0]) + 1im*conj(imagfringe[l,0])
    end
    for m = 1:mmax
        block = blocks[m+1]
        α1 = 2α-1 # positive m
        for l = m:lmax
            block[l-m+1,α1] = conj(realfringe[l,m]) + 1im*conj(imagfringe[l,m])
        end
        α2 = 2α-0 # negative m
        for l = m:lmax
            block[l-m+1,α2] = conj(realfringe[l,m]) - 1im*conj(imagfringe[l,m])
        end
    end
end

function setindex!(transfermatrix::TransferMatrix, block, m, channel)
    ν = transfermatrix.frequencies[channel]
    filename = block_filename(m, ν)
    open(joinpath(transfermatrix.path, filename), "w") do file
        write(file, size(block, 2), size(block, 1), block.')
    end
    block
end

function getindex(transfermatrix::TransferMatrix, m, channel)
    local block
    ν = transfermatrix.frequencies[channel]
    filename = block_filename(m, ν)
    open(joinpath(transfermatrix.path, filename), "r") do file
        sz = tuple(read(file, Int, 2)...)
        block = read(file, Complex128, sz)
    end
    block.'
end

#=
doc"""
    preserve_singular_values(B::TransferMatrix)

Construct a matrix that projects the $m$-modes onto a lower dimensional
space while preserving all the singular values of the transfer matrix.

Multiplying by this matrix will compress the data, make the transfer
matrix square, and leave the information about the sky untouched.
"""
function preserve_singular_values(B::TransferMatrix)
    N = length(B.blocks)
    blocks = Array{MatrixBlock}(N)
    for i = 1:N
        U,σ,V = svd(B.blocks[i])
        blocks[i] = MatrixBlock(U')
    end
    Blocks(blocks)
end
=#

