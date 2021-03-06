# Copyright (c) 2015-2017 Michael Eastwood
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

using Documenter, BPJSpec

makedocs(
    format = :html,
    sitename  = "BPJSpec",
    authors   = "Michael Eastwood",
    linkcheck = ("local" in ARGS),
    html_prettyurls = !("local" in ARGS),
    pages = [
        "Home"                 => "index.md",
        "Noise Model"          => "noise-model.md",
        "Foreground Filtering" => "foreground-filtering.md",
        "Imaging"              => "imaging.md",
        "Power Spectrum Estimation" => "power-spectrum-estimation.md",
        "Enormous Matrices"    => "enormous-matrices.md",
        "Wrappers"             => "wrappers.md"
    ]
)

deploydocs(
    repo   = "github.com/mweastwood/BPJSpec.jl.git",
    julia  = "0.6",
    osname = "linux",
    target = "build",
    deps   = nothing,
    make   = nothing
)

