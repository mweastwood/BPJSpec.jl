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

immutable TransferMeta <: Metadata
    lmax::Int
    m::UnitRange{Int}
    ν::Vector{Float64}

    function TransferMeta(lmax,m,ν)
        if length(m) > 1 && length(ν) > 1
            error("Cannot simultaneously have multiple values of m and multiple frequency channels.")
        end
        new(lmax,m,ν)
    end
end

TransferMeta(lmax::Int,m::Int,ν::AbstractVector) = TransferMeta(lmax,m,collect(ν))
TransferMeta(lmax::Int,mmax::Int,ν::Float64) = TransferMeta(lmax,0:mmax,[ν])

==(lhs::TransferMeta,rhs::TransferMeta) = lhs.lmax == rhs.lmax && lhs.m == rhs.m && lhs.ν == rhs.ν

doc"""
    typealias TransferMatrix Blocks{MatrixBlock, TransferMeta}

This type stores the transfer matrix organized into blocks.

The transfer matrix represents the instrumental response of the
interferometer to the spherical harmonic coefficients of the sky.
This matrix is usually very large, but we can use its block
diagonal structure to work with parts of the matrix separately.
"""
typealias TransferMatrix Blocks{MatrixBlock, TransferMeta}

initial_block_size(::Type{TransferMatrix}, Nbase, lmax, m) = (two(m)*Nbase, lmax-m+1)

function call(::Type{TransferMatrix}, Nbase::Int, lmax::Int, mmax::Int, ν::Float64)
    meta = TransferMeta(lmax,mmax,ν)
    blocks = MatrixBlock[]
    for m = 0:mmax
        sz = initial_block_size(TransferMatrix,Nbase,lmax,m)
        push!(blocks,MatrixBlock(sz))
    end
    TransferMatrix(blocks,meta)
end

is_single_frequency(meta::TransferMeta) = length(meta.ν) == 1
is_single_m(meta::TransferMeta) = length(meta.m) == 1
is_single_frequency(B::TransferMatrix) = is_single_frequency(B.meta)
is_single_m(B::TransferMatrix) = is_single_m(B.meta)

lmax(meta::TransferMeta) = meta.lmax
lmax(B::TransferMatrix) = lmax(B.meta)

mmax(meta::TransferMeta) = maximum(meta.m)
mmax(B::TransferMatrix) = mmax(B.meta)

Nfreq(meta::TransferMeta) = length(meta.ν)
Nfreq(B::TransferMatrix) = length(B.meta)

"""
    transfer(ms::MeasurementSet, beam::TTCal.BeamModel, channel; lmax = 100, mmax = 100)

Generate a transfer matrix where the instrumental parameters are read from a
measurement set.
"""
function transfer(ms::MeasurementSet, beam::TTCal.BeamModel, channel;
                  lmax::Int = 100, mmax::Int = 100)
    healpix_beam = itrf_beam(ms.frame,beam,ms.ν[channel])
    u,v,w = itrf_baselines(ms)
    phasecenter = itrf_phasecenter(ms)
    transfer(healpix_beam, u, v, w, ms.ν[channel],
             phasecenter, lmax=lmax, mmax=mmax)
end

"""
    transfer(beam::HealpixMap, u, v, w, ν, phasecenter; lmax = 100, mmax = 100)

Construct a transfer matrix for the given beam model and observational parameters.
"""
function transfer(beam::HealpixMap, u, v, w, ν, phasecenter;
                  lmax::Int = 100, mmax::Int = 100)
    Nbase = length(u)
    B = TransferMatrix(Nbase,lmax,mmax,ν)
    transfer!(B,beam,u,v,w,ν,phasecenter,lmax=lmax,mmax=mmax)
    B
end

function transfer!(B::TransferMatrix,
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
    l = ReentrantLock()
    increment_progress() = (lock(l); next!(p); unlock(l))
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
        B[1][α,l+1] = conj(realfringe[l,0]) + 1im*conj(imagfringe[l,0])
    end
    for m = 1:mmax, l = m:lmax
        α1 = α         # positive m
        α2 = α + Nbase # negative m
        B[m+1][α1,l-m+1] = conj(realfringe[l,m]) + 1im*conj(imagfringe[l,m])
        B[m+1][α2,l-m+1] = conj(realfringe[l,m]) - 1im*conj(imagfringe[l,m])
    end
end

function save(filename, B::TransferMatrix)
    if !isfile(filename)
        jldopen(filename,"w",compress=true) do file
            file["description"] = "transfer matrix"
        end
    end

    jldopen(filename,"r+",compress=true) do file
        if read(file["description"]) != "transfer matrix"
            error("Attempting to write to a file that does not contain a transfer matrix.")
        end

        if is_single_frequency(B)
            ν = B.meta.ν[1]
            name = @sprintf("%.3fMHz",ν/1e6)
            name in names(file) || g_create(file,name)
            group = file[name]
            for m = 0:mmax(B)
                block = B[m+1]
                group[string(m)] = block.block
            end
        elseif is_single_m(B)
            m = B.meta.m[1]
            for β = 1:Nfreq(B)
                ν = B.meta.ν[β]
                name = @sprintf("%.3fMHz",ν/1e6)
                name in names(file) || g_create(file,name)
                group = file[name]
                block = B[β]
                group[string(m)] = block.block
            end
        end
    end
end

function load(filename, meta::TransferMeta)
    blocks = MatrixBlock[]
    jldopen(filename,"r") do file
        if read(file["description"]) != "transfer matrix"
            error("Attempting to read from a file that does not contain a transfer matrix.")
        end

        if is_single_frequency(meta)
            ν = meta.ν[1]
            name = @sprintf("%.3fMHz",ν/1e6)
            group = file[name]
            for m = 0:mmax(meta)
                block = group[string(m)] |> read
                push!(blocks,MatrixBlock(block))
            end
        elseif is_single_m(meta)
            m = meta.m[1]
            for β = 1:Nfreq(meta)
                ν = meta.ν[β]
                name = @sprintf("%.3fMHz",ν/1e6)
                group = file[name]
                block = group[string(m)] |> read
                push!(blocks,MatrixBlock(block))
            end
        end
    end
    TransferMatrix(blocks,meta)
end

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

