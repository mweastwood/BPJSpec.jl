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

abstract VectorBlock # block that composes a vector
abstract MatrixBlock
typealias Block Union{VectorBlock,MatrixBlock}

Base.size(A::Block) = size(A.block)
Base.length(A::Block) = length(A.block)
Base.svd(A::Block) = svd(A.block)

abstract BlockVector # a vector composed of blocks
abstract BlockDiagonalMatrix

