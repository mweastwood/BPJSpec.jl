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
    itrf_baselines(ms::MeasurementSet)

Calculate the $(u,v,w)$ coordinates in the ITRF coordinate system with the phase
center at true north.
"""
function itrf_baselines(ms::MeasurementSet)
    antenna_table = Table(ms.table[kw"ANTENNA"])
    positions = antenna_table["POSITION"]
    unlock(antenna_table)

    # calculate the baselines while skipping the autocorrelations
    Nbase = ms.Nbase - ms.Nant
    u = zeros(Nbase)
    v = zeros(Nbase)
    w = zeros(Nbase)
    idx = 1
    for α = 1:ms.Nbase
        ms.ant1[α] == ms.ant2[α] && continue
        u[idx] = positions[1,ms.ant1[α]] - positions[1,ms.ant2[α]]
        v[idx] = positions[2,ms.ant1[α]] - positions[2,ms.ant2[α]]
        w[idx] = positions[3,ms.ant1[α]] - positions[3,ms.ant2[α]]
        idx += 1
    end
    u,v,w
end

doc"""
    itrf_phasecenter(ms::MeasurementSet)

Calculate the $(l,m,n)$ coordinates corresponding to the actual phase
center of the array with respsect to true north.
"""
function itrf_phasecenter(ms::MeasurementSet)
    itrf = measure(ms.frame,ms.phase_direction,dir"ITRF")
    vector(itrf)
end

doc"""
    itrf_beam(frame::ReferenceFrame, beam::BeamModel, frequency)

Create a Healpix map of a TTCal beam model.
At the moment only the $I \rightarrow I$ element of the Mueller
matrix is used.
"""
function itrf_beam(frame::ReferenceFrame, beam::TTCal.BeamModel, frequency)
    zenith = Direction(dir"AZEL",Quantity(0.0,"deg"),Quantity(90.0,"deg"))
    zenith_itrf = measure(frame,zenith,dir"ITRF")
    zenith_vec  = [vector(zenith_itrf)...]
    global_north = [0.0,0.0,1.0]
    local_north  = gramschmidt(global_north,zenith_vec)
    map = HealpixMap(zeros(nside2npix(512)))
    @showprogress 1 "Creating map of the beam..." for i = 1:length(map)
        vec = LibHealpix.pix2vec_ring(512,i)
        el  = π/2 - angle_between(vec,zenith_vec)
        vec = gramschmidt(vec,zenith_vec)
        az  = angle_between(vec,local_north)
        if el < 0
            map[i] = 0
        else
            J = beam(frequency,az,el)
            M = MuellerMatrix(J)
            map[i] = M.mat[1,1]
        end
    end
    map
end

doc"""
    immutable TransferMatrixBlock <: MatrixBlock

This type stores a single block of the transfer matrix.

**Fields:**

* `block` stores the actual matrix block
* `lmax` stores the maximum value of $l$
* `m` stores the value of $m$ corresponding to this block of the matrix
* `ν` stores the frequency (in Hz)
"""
immutable TransferMatrixBlock <: MatrixBlock
    block::Matrix{Complex128}
    lmax::Int
    m::Int
    ν::Float64
end

default_size(::Type{TransferMatrixBlock},Nbase,lmax,m) = (two(m)*Nbase, lmax-m+1)

"""
    TransferMatrixBlock(Nbase,lmax,m,ν)

Create an empty transfer matrix block of the default size.
"""
function TransferMatrixBlock(Nbase,lmax,m,ν)
    sz = default_size(TransferMatrixBlock,Nbase,lmax,m)
    block = zeros(Complex128,sz)
    TransferMatrixBlock(block,lmax,m,ν)
end

@enum Representation one_ν one_m

doc"""
    immutable TransferMatrix{rep}

The transfer matrix represents the instrumental response of the
interferometer to the spherical harmonic coefficients of the sky.

The transfer matrix is usually very large, but we can use its block
diagonal structure to work with parts of the matrix separately.
We can choose between:

1. Working with one frequency channel at a time, but all values of $m$,
   where $m$ is the azimuthal quantum number of the spherical harmonics, or
2. Working with several frequency channels, but only one value of $m$.

The former representation is more useful while calculating elements of
the transfer matrix, but the latter representation is more useful for
foreground filtering and cosmological power spectrum estimation (where
frequency structure is important). The type parameter `sym` is therefore
the switch that indicates whether we are working with the first
representation or the second.

    TransferMatrix{one_ν} # one frequency channel, many m
    TransferMatrix{one_m} # many frequency channels, one m

**Fields:**

* `blocks` is a list of `TransferMatrixBlock`s
"""
immutable TransferMatrix{rep}
    blocks::Vector{TransferMatrixBlock}
end

"""
    TransferMatrix(blocks)

Construct a transfer matrix from the list of blocks.

The representation (`one_ν` or `one_m`) is inferred from the provided
list of blocks.
"""
function TransferMatrix(blocks)
    lmax = [block.lmax for block in blocks]
    length(unique(lmax)) == 1 || error("All blocks in a TransferMatrix must have a consistent lmax.")

    m = [block.m for block in blocks]
    ν = [block.ν for block in blocks]
    mmax = length(m) - 1

    # if all the blocks are the same frequency and are sorted in order
    # of increasing m, we have the one_ν representation
    if length(unique(ν)) == 1 && m == collect(0:mmax)
        return TransferMatrix{one_ν}(blocks)
    end

    # if all the blocks have different frequencies (not necessarily sorted)
    # but only have one value of m, we have the one_m
    if length(unique(ν)) == length(ν) && length(unique(m)) == 1
        return TransferMatrix{one_m}(blocks)
    end

    error("""
    The list of blocks must either:

    * have a consistent frequency channel representing all values of m from 0 to some mmax, or
    * have different frequency channels, but a consistent value of m
    """)
end

function TransferMatrix(Nbase,lmax,mmax,ν::Float64)
    [TransferMatrixBlock(Nbase,lmax,m,ν) for m = 0:mmax] |> TransferMatrix
end

# α labels the baseline
# β labels the frequency channel
# l and m label the spherical harmonic

getindex(B::TransferMatrixBlock,α,l) = B.block[α,l-B.m+1]
setindex!(B::TransferMatrixBlock,x,α,l) = B.block[α,l-B.m+1] = x

getindex(B::TransferMatrix{one_ν},m) = B.blocks[m+1]
getindex(B::TransferMatrix{one_ν},α,l,m) = B[m][α,l]
setindex!(B::TransferMatrix{one_ν},x,α,l,m) = B[m][α,l] = x

getindex(B::TransferMatrix{one_m},β) = B.blocks[β]
getindex(B::TransferMatrix{one_m},α,β,l) = B[β][α,l]
setindex!(B::TransferMatrix{one_m},x,α,β,l) = B[β][α,l] = x

frequency(B::TransferMatrix{one_ν}) = B[0].ν
lmax(B::TransferMatrix{one_ν}) = B[0].lmax
mmax(B::TransferMatrix{one_ν}) = length(B.blocks)-1

function ==(lhs::TransferMatrixBlock, rhs::TransferMatrixBlock)
    lhs.lmax == rhs.lmax && lhs.m == rhs.m && lhs.ν == rhs.ν && lhs.block == rhs.block
end

function =={rep}(lhs::TransferMatrix{rep}, rhs::TransferMatrix{rep})
    lhs.blocks == rhs.blocks
end

#=
function SpectralTransferMatrix(Nbase,Nfreq,lmax,mmax,m)
    blocks = [TransferMatrixBlock(Nbase,lmax,mmax,m) for β = 1:Nfreq]
    SpectralTransferMatrix{lmax,mmax}(m,blocks)
end

Nfreq(B::SpectralTransferMatrix) = length(B.blocks)

function Base.size(B::TransferMatrix)
    x = 0; y = 0
    for m = 0:mmax(B)
        sz = size(B[m])
        x += sz[1]
        y += sz[2]
    end
    x,y
end
function Base.size(B::SpectralTransferMatrix)
    x = 0; y = 0
    for β = 1:Nfreq(B)
        sz = size(B[β])
        x += sz[1]
        y += sz[2]
    end
    x,y
end
=#

"""
    gentransfer(ms::MeasurementSet, beam::TTCal.BeamModel, channel; lmax = 100, mmax = 100)

Generate a transfer matrix where the instrumental parameters are read from a
measurement set.
"""
function gentransfer(ms::MeasurementSet, beam::TTCal.BeamModel, channel;
                     lmax::Int = 100, mmax::Int = 100)
    healpix_beam = itrf_beam(ms.frame,beam,ms.ν[channel])
    u,v,w = itrf_baselines(ms)
    phasecenter = itrf_phasecenter(ms)
    gentransfer(healpix_beam, u, v, w, ms.ν[channel],
                phasecenter, lmax=lmax, mmax=mmax)
end

"""
    gentransfer(beam::HealpixMap, u, v, w, ν, phasecenter; lmax = 100, mmax = 100)

Construct a transfer matrix for the given beam model and observational parameters.
"""
function gentransfer(beam::HealpixMap, u, v, w, ν, phasecenter;
                     lmax::Int = 100, mmax::Int = 100)
    Nbase = length(u)
    B = TransferMatrix(Nbase,lmax,mmax,ν)
    gentransfer!(B,beam,u,v,w,ν,phasecenter,lmax=lmax,mmax=mmax)
    B
end

function gentransfer!(B::TransferMatrix{one_ν},
                      beam::HealpixMap, u, v, w, ν, phasecenter;
                      lmax::Int = 100, mmax::Int = 100)
    Nbase = length(u)
    λ = c / ν
    u = u / λ
    v = v / λ
    w = w / λ

    # distribute the workload across all the available workers
    idx = 1
    nextidx() = (myidx = idx; idx += 1; myidx)
    p = Progress(Nbase, 1, "Dreaming...", 50)
    increment_progress() = next!(p)
    @sync for worker in workers()
        @async while true
            α = nextidx()
            α ≤ Nbase || break
            realfringe,imagfringe = remotecall_fetch(worker,fringes,beam,u[α],v[α],w[α],
                                                                    phasecenter,lmax,mmax)
            pack!(B,realfringe,imagfringe,α,Nbase,lmax,mmax)
            increment_progress()
        end
    end
    B
end

"""
    fringes(beam, u, v, w, phasecenter, lmax, mmax)

Generate the spherical harmonic expansion of the fringe pattern on the sky.
Note that because the Healpix library assumes you are asking for the coefficients
of a real field, there must be one set of coefficients for the real part of
the fringe pattern and one set of coefficients for the imaginary part of the
fringe pattern.
"""
function fringes(beam,u,v,w,phasecenter,lmax,mmax)
    # Account for the extra phase due to the fact that the phase
    # center is not perfectly orthogonal to the baseline.
    extraphase = -2π*(u*phasecenter[1]+v*phasecenter[2]+w*phasecenter[3])
    realfringe,imagfringe = planewave(u,v,w,extraphase,lmax=lmax,mmax=mmax)
    # Use spherical harmonic transforms to incorporate the effects of the beam.
    realfringe = map2alm(beam.*alm2map(realfringe,nside=512),lmax=lmax,mmax=mmax)
    imagfringe = map2alm(beam.*alm2map(imagfringe,nside=512),lmax=lmax,mmax=mmax)
    realfringe,imagfringe
end

"""
    pack!(B, realfringe, imagfringe, α, Nbase, lmax, mmax)

Having calculated the spherical harmonic expansion of the fringe pattern,
pack those numbers into the transfer matrix.
"""
function pack!(B,realfringe,imagfringe,α,Nbase,lmax,mmax)
    # Note that all the conjugations in this function come about because
    # Shaw et al. 2014, 2015 expand the fringe pattern in terms of the
    # spherical harmonic conjugates while we've expanded the fringe pattern
    # in terms of the spherical harmonics.
    for l = 0:lmax
        B[α,l,0] = conj(realfringe[l,0]) + 1im*conj(imagfringe[l,0])
    end
    for m = 1:mmax, l = m:lmax
        α1 = α         # positive m
        α2 = α + Nbase # negative m
        B[α1,l,m] = conj(realfringe[l,m]) + 1im*conj(imagfringe[l,m])
        B[α2,l,m] = conj(realfringe[l,m]) - 1im*conj(imagfringe[l,m])
    end
end

function save_transfermatrix(filename, B::TransferMatrix{one_ν})
    if !isfile(filename)
        jldopen(filename,"w",compress=true) do file
            file["lmax"] = lmax(B)
            file["mmax"] = mmax(B)
        end
    end

    jldopen(filename,"r+",compress=true) do file
        read(file["lmax"]) == lmax(B) || error("lmax is inconsistent")
        read(file["mmax"]) == mmax(B) || error("mmax is inconsistent")

        name_ν  = @sprintf("%.3fMHz",frequency(B)/1e6)
        name_ν in names(file) && error("group $(name_ν) already exists in file")
        group_ν = g_create(file,name_ν)

        for m = 0:mmax(B)
            group_m = g_create(group_ν,string(m))
            group_m["block"] = B[m].block
        end
    end
end

function load_transfermatrix(filename, frequency)
    blocks = TransferMatrixBlock[]
    jldopen(filename,"r") do file
        lmax = file["lmax"] |> read
        mmax = file["mmax"] |> read

        name_ν  = @sprintf("%.3fMHz",frequency/1e6)
        group_ν = file[name_ν]

        for m = 0:mmax
            group_m = group_ν[string(m)]
            block = group_m["block"] |> read
            push!(blocks,TransferMatrixBlock(block,lmax,m,frequency))
        end
    end
    TransferMatrix(blocks)
end

#=
"""
    SpectralTransferMatrix(m,matrices::Vector{TransferMatrix})

Construct a SpectralTransferMatrix from the given list of transfer matrices.
Each transfer matrix should correspond to a different frequency channel.
"""
function SpectralTransferMatrix(m,matrices::Vector{TransferMatrix})
    lmax′ = lmax(matrices[1])
    mmax′ = mmax(matrices[1])
    blocks = TransferMatrixBlock[]
    for matrix in matrices
        if lmax(matrix) != lmax′ || mmax(matrix) != mmax′
            error("The transfer matrices must all have the same lmax and mmax.")
        end
        push!(blocks,matrix[m])
    end
    SpectralTransferMatrix{lmax′,mmax′}(m,blocks)
end

function Base.full(B::SpectralTransferMatrix)
    out = zeros(Complex128,size(B))
    idx1 = 1; idx2 = 1
    for β = 1:Nfreq(B)
        block = B[β].block
        out[idx1:idx1+size(block,1)-1,
            idx2:idx2+size(block,2)-1] = block
        idx1 += size(block,1)
        idx2 += size(block,2)
    end
    out
end
=#

