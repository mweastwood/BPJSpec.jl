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

struct AngularCovarianceMetadata <: MatrixMetadata
    lmax        :: Int
    frequencies :: Vector{typeof(1.0*u"Hz")}
    bandwidth   :: Vector{typeof(1.0*u"Hz")}
end
function Base.show(io::IO, metadata::AngularCovarianceMetadata)
    @printf(io, "{lmax: %d, ν: %.3f MHz..%.3f MHz, Δν: %.3f MHz total}",
            metadata.lmax,
            ustrip(uconvert(u"MHz", metadata.frequencies[1])),
            ustrip(uconvert(u"MHz", metadata.frequencies[end])),
            ustrip(uconvert(u"MHz", sum(metadata.bandwidth))))
end
indices(metadata::AngularCovarianceMetadata) = L(0):L(metadata.lmax)
number_of_blocks(metadata::AngularCovarianceMetadata) = metadata.lmax+1

"""
Generally speaking we want to index all of our matrices by the integer m. However, our angular
covariance matrices are block-diagonal in l, which is also an integer. We will therefore use this
type to indicate that we'd like to index with l and indexing with an integer will continue to mean
indexing with m.
"""
struct L <: Integer
    l :: Int
end
Base.convert(::Type{L}, l::Int) = L(l)
Base.convert(::Type{Int}, l::L) = l.l
Base.promote_rule(::Type{L}, ::Type{Int}) = L
Base.oneunit(::L) = L(1)
Base.:≤(lhs::L, rhs::L) = lhs.l ≤ rhs.l
Base.:+(lhs::L, rhs::L) = L(lhs.l + rhs.l)
Base.:-(lhs::L, rhs::L) = L(lhs.l - rhs.l)

struct AngularCovarianceMatrix{S} <: AbstractBlockMatrix
    matrix :: BlockMatrix{Matrix{Float64}, 1, AngularCovarianceMetadata, S}
end
function AngularCovarianceMatrix(storage::Mechanism, mmax, frequencies, bandwidth)
    metadata = AngularCovarianceMetadata(mmax, frequencies, bandwidth)
    matrix = BlockMatrix{Matrix{Float64}, 1}(storage, metadata)
    AngularCovarianceMatrix(matrix)
end
function AngularCovarianceMatrix(path::String)
    AngularCovarianceMatrix(BlockMatrix{Matrix{Float64}, 1}(path))
end

Base.getindex(matrix::AngularCovarianceMatrix, l::L) = matrix.matrix[l+1]
Base.getindex(matrix::AngularCovarianceMatrix, l, m) = matrix.matrix[l+1]
Base.setindex!(matrix::AngularCovarianceMatrix, block, l) = matrix.matrix[l+1] = block

function Base.getindex(matrix::AngularCovarianceMatrix, m)
    blocks = [matrix[l, m] for l = m:lmax(matrix)]
    Nfreq  = length(frequencies(matrix))
    Nl     = lmax(matrix)-m+1
    Ntotal = Nfreq*Nl
    output = zeros(Float64, Ntotal, Ntotal)
    for l = m:lmax(matrix)
        block = blocks[l-m+1]
        for β1 = 1:Nfreq, β2 = 1:Nfreq
            idx1 = (β1-1)*Nl + l - m + 1
            idx2 = (β2-1)*Nl + l - m + 1
            output[idx1, idx2] = block[β1, β2]
        end
    end
    output
end

function compute!(matrix::AngularCovarianceMatrix, component::SkyComponent; progress=false)
    ν  = frequencies(matrix)
    Δν =   bandwidth(matrix)
    Nfreq = length(ν)
    args  = precomputation(component)
    progress && (prg = Progress(lmax(matrix)+1))
    for l = 0:lmax(matrix)
        block = zeros(Float64, Nfreq, Nfreq)
        for β1 = 1:Nfreq
            ν1  =  ν[β1]
            Δν1 = Δν[β1]
            block[β1, β1] = compute(component, l, ν1, Δν1, ν1, Δν1, args)
            for β2 = β1+1:Nfreq
                ν2  =  ν[β2]
                Δν2 = Δν[β2]
                block[β1, β2] = compute(component, l, ν1, Δν1, ν2, Δν2, args)
                block[β2, β1] = block[β1, β2]
            end
        end
        matrix[l] = block
        progress && next!(prg)
    end
end

precomputation(component) = nothing

"Compute the covariance while accounting for bandwidth smearing."
function compute(component, l, ν1, Δν1, ν2, Δν2, args = nothing)
    if args == nothing
        integrand = x -> ustrip(uconvert(u"K^2", component(l, x[1]*u"Hz", x[2]*u"Hz")))
    else
        integrand = x -> ustrip(uconvert(u"K^2", component(l, x[1]*u"Hz", x[2]*u"Hz", args)))
    end
    xmin = ustrip.(uconvert.(u"Hz", [ν1-Δν1/2, ν2-Δν2/2]))
    xmax = ustrip.(uconvert.(u"Hz", [ν1+Δν1/2, ν2+Δν2/2]))

    # We choose pcubature because usually the functions are pretty smooth and pcubature tends to
    # converge 100x faster than hcubature in this usual case.  However, We need to limit the number
    # of evaluations here otherwise the adaptive routine tends to get stuck trying to achieve the
    # tolerance we asked for. This tends to happen when the frequency smearing straddles the border
    # between two kperp cells (due to the frequency dependent distance).
    #
    # Recall: kperp = l / χ and χ is a function of the redshift (equivalently frequency).

    val, err = pcubature(integrand, xmin, xmax, reltol=1e-8, maxevals=100_000)
    val / ustrip(uconvert(u"Hz^2", Δν1*Δν2))
end

