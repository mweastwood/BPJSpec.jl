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

# This function is useful to handle some of the special casing required for m == 0
two(m) = ifelse(m != 0, 2, 1)

# These functions define the filenames used when blocks are written to disk
block_filename(ν) = @sprintf("%.6fMHz.block", ν/1e6)
block_filename(m, ν) = @sprintf("%.6fMHz-%04d.block", ν/1e6, m)

# This function defines the ordering for blocks of a given m and frequency channel
# This ordering is assumed at various points so a careful check is needed before changing
block_index(mmax, m, channel) = (mmax+1)*(channel-1) + m + 1

"""
    @define(expr)

We use this macro to automatically generate some helper methods for
the type definitions below.
"""
macro define(expr)
    type_name = expr.args[2]
    field_names = Symbol[]
    for field in expr.args[3].args
        if field.head == :(::)
            push!(field_names, field.args[1])
        elseif field.head == :line
            # ignore line info
        elseif field.head == :function
            # ignore inner constructors
        else
            error("check type definition ast parsing ($(field.head))")
        end
    end

    isequal_expr = :(true)
    for field in field_names
        isequal_expr = :($isequal_expr && lhs.$field == rhs.$field)
    end

    quote
        Base.@__doc__ $expr
        function ==(lhs::$type_name, rhs::$type_name)
            $isequal_expr
        end
    end |> esc
end

doc"""
    GriddedVisibilities

This type represents visibilities on a sidereal time grid.

# Fields

* `path` points to the directory where all of the visibilities are stored
* `Nbase` is the total number of baselines
* `Ntime` is the number of sidereal time grid points
* `frequencies` is the list of frequencies in units of Hz
* `origin` is the sidereal time of the first grid point
* `data` is a list of mmapped arrays where the visibilities are stored
* `weights` is a list of mmapped arrays where the weights are stored

Note that in each of `data` and `weights` there is one mmapped array for each
frequency channel, and each array has dimensions of `Nbase` by `Ntime`.

# Implementation

Note that generally an interferometer integrates for some time $t$
that does not evenly divide a sidereal day. We must therefore pick
some gridding kernel that determines what we do with an integration
that falls between grid points. The gridding kernel currently used
is the triangular function (or hat function), which divides the
visibility amongst the two nearest grid points proportional to its
distance from each grid point.

At the OVRO LWA we picked a 13 second integration time which comes
within one tenth of one second to evenly dividing a sidereal day.
However we cannot guarantee that the correlator starts at a given
sidereal time. Therefore we need to adjust the grid to align with
the integrations. The `origin` parameter is used to accomplish this.

When gridding the visibilities we are going to be making a lot of small writes
to several arrays that may not fit into the system memory. This is why we mmap
the arrays onto the disk.

Experiements suggest that two arrays cannot be mmapped to the same file.
That is instead of being written one after another the two arrays are written
on top of each other. This is why `data` and `weights` are mmapped to two
separate files.
"""
@define immutable GriddedVisibilities
    path :: ASCIIString
    Nbase :: Int
    Ntime :: Int
    frequencies :: Vector{Float64}
    origin :: Float64
    data :: Vector{Matrix{Complex128}}
    weights :: Vector{Matrix{Float64}}
    function GriddedVisibilities(path, Nbase, Ntime, frequencies, origin, data, weights)
        0 ≤ origin < 1 || throw(ArgumentEttor("The sidereal time must be in the interval [0,1)"))
        new(path, Nbase, Ntime, frequencies, origin, data, weights)
    end
end

doc"""
    MModes

This type represents the $m$-modes measured by an interferometer.

An $m$-mode is the Fourier transform of a visibility with
respect to sidereal time. These are related to the spherical
harmonic coefficients of the sky brightness through a matrix equation.

```math
v = Ba,
```

where $v$ is the vector of $m$-modes, $B$ is the transfer
matrix, and $a$ is the vector of spherical harmonic coefficients.

# Fields

* `path` points to the directory that contains the blocks of $m$-modes
* `mmax` is the maximum value of the azimuthal quantum number $m$
* `frequencies` is the list of frequencies in units of Hz
* `blocks` is a list of mmapped vectors where the actual $m$-modes are stored

Note that there is one mmapped vector for each frequency and each value of $m$.
"""
@define immutable MModes
    path :: ASCIIString
    mmax :: Int
    frequencies :: Vector{Float64}
    blocks :: Vector{Vector{Complex128}}
end

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
* `frequencies` is the list of frequencies in units of Hz

# Implementation

Because the transfer matrix can be incredibly large it is almost
certainly the case that it cannot fit into the system memory.
However each individual block of the matrix should be able to fit.
Therefore the matrix must be stored on disk.
"""
@define immutable TransferMatrix
    path :: ASCIIString
    lmax :: Int
    mmax :: Int
    frequencies :: Vector{Float64}
    function TransferMatrix(path, lmax, mmax, frequencies)
        mmax ≤ lmax || throw(ArgumentError("Transfer matrices require mmax ≤ lmax"))
        new(path, lmax, mmax, frequencies)
    end
end

doc"""
    MultiFrequencyAlm

This type represents a set of spherical harmonic coefficients
for multiple frequency channels.

# Fields

* `path` points to the directory where the spherical harmonic coefficients are stored
* `lmax` is the maximum value of the total angular momentum quantum number $l$
* `mmax` is the maximum value of the azimuthal quantum number $m$
* `frequencies` is the list of frequencies in units of Hz
* `alms` is a list of spherical harmonic coefficients $a_{lm}$ for each frequency channel
"""
@define immutable MultiFrequencyAlm
    path :: ASCIIString
    lmax :: Int
    mmax :: Int
    frequencies :: Vector{Float64}
    alms :: Vector{Vector{Complex128}}
end

