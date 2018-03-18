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

struct MultiFrequencyAlm end

function create(::Type{MultiFrequencyAlm}, input::MFBlockVector)
    lmax = mmax = input.mmax
    Nfreq = length(input.frequencies)
    output = create(LMBlockVector, lmax, mmax, input.frequencies, input.bandwidth)
    for m = 0:mmax, l = m:lmax
        output_block = zeros(Complex128, Nfreq)
        for β = 1:Nfreq
            input_block = input[m, β]
            output_block[β] = input_block[l-m+1]
        end
        output[l, m] = output_block
    end
    output
end

function create(::Type{MultiFrequencyAlm}, input::MBlockVector, frequencies, bandwidth)
    lmax = mmax = input.mmax
    Nfreq = length(frequencies)
    output = create(LMBlockVector, lmax, mmax, frequencies, bandwidth)
    for m = 0:mmax, l = m:lmax
        output_block = zeros(Complex128, Nfreq)
        for β = 1:Nfreq
            input_block = input[m]
            output_block[β] = input_block[(lmax-m+1)*(β-1) + (l-m+1)]
        end
        output[l, m] = output_block
    end
    output
end

