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

doc"""
    TransferMatrix

This type represents the transfer matrix of an interferometer.

The transfer matrix represents the instrumental response of the
interferometer to the spherical harmonic coefficients of the sky.
This matrix is usually very large, but we can use its block
diagonal structure to work with parts of the matrix separately.

# Fields

* `path` points to the directory that contains all of the transfer matrix blocks
* `lmax` is the maximum value of the total angular momentum quantum number $l$
* `mmax` is the maximum value of the azimuthal quantum number $m$
* `ν` is the list of frequencies in units of Hz

# Implementation

Because the transfer matrix can be incredibly large it is almost
certainly the case that it cannot fit into the system memory.
However each individual block of the matrix should be able to fit.
Therefore the matrix must be stored on disk.
"""
immutable TransferMatrix
    path :: ASCIIString
    lmax :: Int
    mmax :: Int
    ν :: Vector{Float64}
    function TransferMatrix(path, lmax, mmax, ν)
        mmax ≤ lmax || throw(ArgumentError("Transfer matrices require mmax ≤ lmax"))
        new(path, lmax, mmax, ν)
    end
end

function TransferMatrix(path)
    local lmax, mmax, ν
    open(joinpath(path, "METADATA"), "r") do file
        lmax = read(file, Int)
        mmax = read(file, Int)
        len  = read(file, Int)
        ν    = read(file, Float64, len)
    end
    TransferMatrix(path, lmax, mmax, ν)
end

Nfreq(matrix::TransferMatrix) = length(matrix.ν)
channels(matrix::TransferMatrix) = 1:Nfreq(matrix)

"""
    TransferMatrix(path, meta, lmax, mmax)

Generate the transfer matrix. This will take some time.
"""
function TransferMatrix(path, meta::Metadata, lmax, mmax)
    transfermatrix = TransferMatrix(path, lmax, mmax, meta.channels)
    initialize!(transfermatrix, Nbase(meta))
    generate_transfermatrix!(transfermatrix, meta)
    transfermatrix
end

"""
    initialize!(transfermatrix, Nbase)

Initialize the transfer matrix blocks to be all zeros.
The number of baselines must be specified so that all the blocks have the right size.
"""
function initialize!(transfermatrix::TransferMatrix, Nbase)
    path = transfermatrix.path
    lmax = transfermatrix.lmax
    mmax = transfermatrix.mmax
    ν = transfermatrix.ν
    # create the directory if it doesn't already exist
    isdir(path) || mkdir(path)
    # create the METADATA file to store lmax, mmax, and the list of frequency channels
    open(joinpath(transfermatrix.path, "METADATA"), "w") do file
        write(file, lmax, mmax, length(ν), ν)
    end
end

function generate_transfermatrix!(transfermatrix, meta)
    for ν in transfermatrix.ν
        generate_transfermatrix_onechannel!(transfermatrix, meta, ν)
    end
end

function generate_transfermatrix_onechannel!(transfermatrix, meta, ν)
    lmax = transfermatrix.lmax
    mmax = transfermatrix.mmax
    beam = beam_map(meta, ν)
    frame = TTCal.reference_frame(meta)
    phase_center = measure(frame, meta.phase_center, dir"ITRF")
    # Memory map all the blocks on the master process to avoid having to
    # open/close the files multiple times and to avoid having to read the
    # entire matrix at once.
    blocks = Matrix{Complex128}[]
    for m = 0:transfermatrix.mmax
        filename = block_filename(m, ν)
        open(joinpath(transfermatrix.path, filename), "w+") do file
            sz = (two(m)*Nbase(meta), lmax-m+1)
            write(file, sz[1], sz[2])
            block = Mmap.mmap(file, Matrix{Complex128}, sz)
            push!(blocks, block)
        end
    end
    @distribute for α = 1:Nbase(meta)
        realfringe, imagfringe = @remote fringes(meta, beam, phase_center, lmax, mmax, ν, α)
        pack!(blocks, realfringe, imagfringe, lmax, mmax, Nbase(meta), α)
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
    progress = Progress(length(map), "Creating a Healpix map of the beam: ")
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
        next!(progress)
    end
    map
end

"""
    planewave(u, v, w, phase, lmax, mmax)

Compute the spherical harmonic coefficients corresponding to the
plane wave:

    exp(2im*π*(u*x+v*y+w*z)) * exp(1im*phase)
"""
function planewave(u, v, w, phase, lmax, mmax)
    b = sqrt(u^2 + v^2 + w^2)
    θ = acos(w / b)
    ϕ = atan2(v, u)
    realpart = Alm(Complex128, lmax, mmax)
    imagpart = Alm(Complex128, lmax, mmax)
    for m = 0:mmax, l = m:lmax
        alm1 = 4π * (1im)^l * j(l,2π*b) * conj(       Y(l,+m,θ,ϕ)) * exp(1im*phase)
        alm2 = 4π * (1im)^l * j(l,2π*b) * conj((-1)^m*Y(l,-m,θ,ϕ)) * exp(1im*phase)
        realpart[l,m] = (alm1 + conj(alm2))/2
        imagpart[l,m] = (alm1 - conj(alm2))/2im
    end
    realpart, imagpart
end

"""
    fringes(meta, beam, phase_center, lmax, mmax, ν, α)

Generate the spherical harmonic expansion of the fringe pattern on the sky.

Note that because the Healpix library assumes you are asking for the coefficients
of a real field, there must be one set of coefficients for the real part of
the fringe pattern and one set of coefficients for the imaginary part of the
fringe pattern.
"""
function fringes(meta, beam, phase_center, lmax, mmax, ν, α)
    antenna1 = meta.antennas[meta.baselines[α].antenna1]
    antenna2 = meta.antennas[meta.baselines[α].antenna2]
    λ = c / ν
    u = (antenna1.position.x - antenna2.position.x) / λ
    v = (antenna1.position.y - antenna2.position.y) / λ
    w = (antenna1.position.z - antenna2.position.z) / λ
    if hypot(hypot(u, v), w) < eps(Float64)
        # Special case for the autocorrelations to prevent NaNs from entering into
        # the transfer matrix elements.
        realfringe = map2alm(beam, lmax, mmax, iterations=0)
        imagfringe = Alm(Complex128, lmax, mmax) # all zeros
    else
        # Account for the extra phase due to the fact that the phase center is
        # not perfectly orthogonal to the baseline.
        extraphase = -2π*(u*phase_center.x + v*phase_center.y + w*phase_center.z)
        realfringe, imagfringe = planewave(u, v, w, extraphase, lmax, mmax)
        # Use spherical harmonic transforms to incorporate the effects of the beam.
        realfringe = map2alm(beam .* alm2map(realfringe, 512), lmax, mmax, iterations=0)
        imagfringe = map2alm(beam .* alm2map(imagfringe, 512), lmax, mmax, iterations=0)
    end
    realfringe, imagfringe
end

"""
    pack!(blocks, realfringe, imagfringe, lmax, mmax, Nbase, α)

Having calculated the spherical harmonic expansion of the fringe pattern,
pack those numbers into the transfer matrix.
"""
function pack!(blocks, realfringe, imagfringe, lmax, mmax, Nbase, α)
    # Note that all the conjugations in this function come about because
    # Shaw et al. 2014, 2015 expand the fringe pattern in terms of the
    # spherical harmonic conjugates while we've expanded the fringe pattern
    # in terms of the spherical harmonics.
    for l = 0:lmax
        blocks[1][α,l+1] = conj(realfringe[l,0]) + 1im*conj(imagfringe[l,0])
    end
    for m = 1:mmax, l = m:lmax
        α1 = α         # positive m
        α2 = α + Nbase # negative m
        blocks[m+1][α1,l-m+1] = conj(realfringe[l,m]) + 1im*conj(imagfringe[l,m])
        blocks[m+1][α2,l-m+1] = conj(realfringe[l,m]) - 1im*conj(imagfringe[l,m])
    end
end

function setindex!(transfermatrix, block, m, channel)
    ν = transfermatrix.ν[channel]
    filename = block_filename(m, ν)
    open(joinpath(transfermatrix.path, filename), "w") do file
        write(file, size(block, 1), size(block, 2), block)
    end
    block
end

function getindex(transfermatrix::TransferMatrix, m, channel)
    local block
    ν = transfermatrix.ν[channel]
    filename = block_filename(m, ν)
    open(joinpath(transfermatrix.path, filename), "r") do file
        sz = tuple(read(file, Int, 2)...)
        block = read(file, Complex128, sz)
    end
    block
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

