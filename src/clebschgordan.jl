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

if !isfile(joinpath(dirname(@__FILE__),"cg-coeff.csv"))
    error("""
          Could not find "cg-coeff.csv".

          BPJSpec relies on a table of Clebsch-Gordan coefficients that
          can be generated with Mathematica.

          Please run the following command in Mathematica and move the
          generated "cg-coeff.csv" file to BPJSpec's src/ directory.

              Export["cg-coeff.csv", 
                     Flatten[Table[N[
                     ClebschGordan[{l1, m1}, {l2, m2}, {L, m1 + m2}]],
                     {l1, 0, 10}, {m1, -l1, l1},
                     {l2, 0, 100}, {m2, -l2, l2},
                     {L, Max[Abs[l1 - l2], Abs[m1 + m2]], l1 + l2}]]]

          Be aware that this command will take several hours to complete.
          However it only needs to be run once.
          """)
end

immutable CGTable
    coeff::Vector{Float64} # one massive flattened list of all the CG coefficients
    zerom::Matrix{Vector{Float64}} # all the CG coefficients with m1 == m2 == 0
                                   # These coefficients are stored separately to
                                   # make finding them easier (I don't know how to
                                   # index into the flattened list).
end

# Do not change these parameters without updating the error message above.
l1max(::Type{CGTable}) = 10
l2max(::Type{CGTable}) = 100

function CGTable()
    coeff = squeeze(readcsv(joinpath(dirname(@__FILE__),"cg-coeff.csv")),2)

    # Find the coefficients with m1 == m2 == 0
    zerom = [Float64[] for l1 = 0:10, l2 = 0:100]
    idx = 1
    for l1 = 0:l1max(CGTable), m1 = -l1:l1,
        l2 = 0:l2max(CGTable), m2 = -l2:l2,
        L = max(abs(l1-l2),abs(m1+m2)):l1+l2

        if m1 == m2 == 0
            push!(zerom[l1+1,l2+1],coeff[idx])
        end
        idx += 1
    end

    CGTable(coeff,zerom)
end

@doc """
Multiply two spherical harmonic expansions together
by means of the Clebsch-Gordan series.

http://functions.wolfram.com/Polynomials/SphericalHarmonicY/16/ShowAll.html
""" ->
function multiply(alm1::Alm, alm2::Alm, cg::CGTable)
    output = Alm(Complex128,max(lmax(alm1),lmax(alm2)),
                            max(mmax(alm1),mmax(alm2)))
    multiply!(output,alm1,alm2,cg)
    output
end

function multiply!(output::Alm, alm1::Alm, alm2::Alm, cg::CGTable)
    idx = 1
    for l1 = 0:l1max(CGTable), m1 = -l1:l1,
        l2 = 0:l2max(CGTable), m2 = -l2:l2,
        L = max(abs(l1-l2),abs(m1+m2)):l1+l2

        M = m1 + m2
        if (L ≤ lmax(output) && 0 ≤ M ≤ mmax(output)
            && l1 ≤ lmax(alm1) && m1 ≤ mmax(alm1)
            && l2 ≤ lmax(alm2) && m2 ≤ mmax(alm2))

            idx2 = L - abs(l1-l2) + 1

            output[L,abs(M)] += (sqrt((2l1+1)*(2l2+1)/(4π*(2L+1)))
                                    * cg.coeff[idx]
                                    * cg.zerom[l1+1,l2+1][idx2]
                                    * alm1[l1,m1]
                                    * alm2[l2,m2])
        end
        idx += 1
    end
end

