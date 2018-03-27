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

struct AngularCovarianceMatrix end

function create(::Type{AngularCovarianceMatrix}, storage::Mechanism, component::SkyComponent,
                lmax, frequencies, bandwidth;
                rm=false, progress=false)
    output = create(LBlockMatrix, storage, lmax, frequencies, bandwidth, rm=rm)
    compute!(AngularCovarianceMatrix, output, component, progress=progress)
    output
end

function compute!(::Type{AngularCovarianceMatrix}, matrix::LBlockMatrix,
                  component::SkyComponent; progress=false)
    ν  = matrix.frequencies
    Δν = matrix.bandwidth
    Nfreq = length(ν)
    args  = precomputation(component)
    progress && (prg = Progress(matrix.lmax+1))
    for l = L(0):L(matrix.lmax)
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

