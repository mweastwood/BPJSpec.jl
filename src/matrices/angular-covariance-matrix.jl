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

struct AngularCovarianceMatrix <: BlockMatrix
    path        :: String
    progressbar :: Bool
    distribute  :: Bool
    cached      :: Ref{Bool}
    lmax        :: Int
    frequencies :: Vector{typeof(1.0*u"Hz")}
    bandwidth   :: Vector{typeof(1.0*u"Hz")}
    component   :: SkyComponent
    blocks      :: Vector{Matrix{Float64}}

    function AngularCovarianceMatrix(path, lmax, frequencies, bandwidth, component, write=true;
                                     progressbar=false, distribute=false, cached=false)
        if write
            isdir(path) || mkpath(path)
            save(joinpath(path, "METADATA.jld2"), "lmax", lmax,
                 "frequencies", frequencies, "bandwidth", bandwidth, "component", component)
        end
        blocks = Matrix{Float64}[]
        output = new(path, progressbar, distribute, Ref(cached),
                     lmax, frequencies, bandwidth, component, blocks)
        if write
            compute!(output)
        elseif cached
            cache!(output)
        end
        output
    end
end

function AngularCovarianceMatrix(path; kwargs...)
    lmax, frequencies, bandwidth, component = load(joinpath(path, "METADATA.jld2"),
                                                   "lmax", "frequencies", "bandwidth", "component")
    AngularCovarianceMatrix(path, lmax, frequencies, bandwidth, component, false; kwargs...)
end

function compute!(matrix::AngularCovarianceMatrix)
    Nfreq = length(matrix.frequencies)
    if progressbar(matrix)
        prg = Progress(matrix.lmax+1)
    end
    args = precomputation(matrix.component)
    for l = 0:matrix.lmax
        block = zeros(Float64, Nfreq, Nfreq)
        for β1 = 1:Nfreq
            ν1  = matrix.frequencies[β1]
            Δν1 = matrix.bandwidth[β1]
            block[β1, β1] = compute(matrix.component, l, ν1, Δν1, ν1, Δν1, args)
            for β2 = β1+1:Nfreq
                ν2 = matrix.frequencies[β2]
                Δν2 = matrix.bandwidth[β2]
                block[β1, β2] = compute(matrix.component, l, ν1, Δν1, ν2, Δν2, args)
                block[β2, β1] = block[β1, β2]
            end
        end
        matrix[l, 0] = block
        progressbar(matrix) && next!(prg)
    end
end

precomputation(component) = nothing

"Compute the covariance while accounting for bandwidth smearing."
function compute(component, l, ν1, Δν1, ν2, Δν2, args = nothing)
    if args == nothing
        integrand = x -> component(l, x[1]*u"Hz", x[2]*u"Hz")
    else
        integrand = x -> component(l, x[1]*u"Hz", x[2]*u"Hz", args)
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

Base.show(io::IO, matrix::AngularCovarianceMatrix) =
    print(io, "AngularCovarianceMatrix: ", matrix.path)

indices(matrix::AngularCovarianceMatrix) =
    [(l, m) for m = 0:matrix.lmax for l = m:matrix.lmax]

function Base.getindex(matrix::AngularCovarianceMatrix, l, m)
    if matrix.cached[]
        return matrix.blocks[l+1]
    else
        return read_from_disk(matrix, l)
    end
end

function Base.setindex!(matrix::AngularCovarianceMatrix, block, l, m)
    if matrix.cached[]
        matrix.blocks[l+1] = block
    else
        write_to_disk(matrix, block, l)
    end
    block
end

function Base.getindex(matrix::AngularCovarianceMatrix, m)
    blocks = [matrix[l, m] for l = m:matrix.lmax]
    Nfreq  = length(matrix.frequencies)
    Nl     = matrix.lmax-m+1
    Ntotal = Nfreq*Nl
    output = zeros(Float64, Ntotal, Ntotal)
    for l = m:matrix.lmax
        block = blocks[l-m+1]
        for β1 = 1:Nfreq, β2 = 1:Nfreq
            idx1 = (β1-1)*Nl + l - m + 1
            idx2 = (β2-1)*Nl + l - m + 1
            output[idx1, idx2] = block[β1, β2]
        end
    end
    output
end

function cache!(matrix::AngularCovarianceMatrix)
    matrix.cached[] = true
    empty!(matrix.blocks)
    for l = 0:matrix.lmax
        push!(matrix.blocks, read_from_disk(matrix, l))
    end
    matrix
end

function flush!(matrix::AngularCovarianceMatrix)
    for l = 0:matrix.lmax
        write_to_disk(matrix, matrix.blocks[l+1], l)
    end
    empty!(matrix.blocks)
    matrix.cached[] = false
    matrix
end

function read_from_disk(matrix::AngularCovarianceMatrix, l::Integer)
    filename   = @sprintf("l=%04d.jld2", l)
    objectname = "block"
    load(joinpath(matrix.path, filename), objectname) :: Matrix{Float64}
end

function write_to_disk(matrix::AngularCovarianceMatrix, block::Matrix{Float64}, l::Integer)
    filename   = @sprintf("l=%04d.jld2", l)
    objectname = "block"
    save(joinpath(matrix.path, filename), objectname, block)
    block
end

