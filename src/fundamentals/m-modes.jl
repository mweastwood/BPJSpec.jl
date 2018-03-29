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
    struct MModes

This singleton type represents the $m$-modes measured by the interferometer. $m$-modes are computed
from a Fourier transform of the measured visibilities over sidereal time $ϕ ∈ [0, 2π)$.

```math
\text{m-mode} = \int_0^{2π} (\text{visibility}) \times e^{-imϕ} \, {\rm d}ϕ
```
"""
struct MModes end

function create(::Type{MModes}, storage::Mechanism,
                metadata::Metadata, hierarchy::Hierarchy; rm=false)
    output = create(MFBlockVector, storage, hierarchy.lmax,
                    metadata.frequencies, metadata.bandwidth, rm=rm)
    output
end

"Compute m-modes from two dimensional matrix of visibilities (time × baseline)."
function compute!(::Type{MModes}, mmodes::MFBlockVector, hierarchy::Hierarchy,
                  visibilities::Matrix{Complex128}, β; dϕ=0.0)
    store!(mmodes, hierarchy, fourier_transform(visibilities), β, dϕ)
end

function fourier_transform(matrix)
    Ntime, Nbase = size(matrix)
    planned_fft = plan_fft(matrix, 1)
    (planned_fft*matrix) ./ Ntime
end

function store!(mmodes, hierarchy, transformed_visibilities, β, dϕ)
    Ntime = size(transformed_visibilities, 1)

    # m = 0
    block = zeros(Complex128, Nbase(hierarchy, 0))
    for (α, α′) in enumerate(baseline_permutation(hierarchy, 0))
        block[α] = transformed_visibilities[1, α′]
    end
    mmodes[0, β] = block

    # m > 0
    for m = 1:mmodes.mmax
        block = zeros(Complex128, 2Nbase(hierarchy, m))
        rotation = cis(m*dϕ)
        for (α, α′) in enumerate(baseline_permutation(hierarchy, m))
            α1 = 2α-1 # positive m
            α2 = 2α-0 # negative m
            block[α1] =      transformed_visibilities[      m+1, α′]  * rotation
            block[α2] = conj(transformed_visibilities[Ntime+1-m, α′]) * rotation
        end
        mmodes[m, β] = block
    end
end

