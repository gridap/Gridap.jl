module CellArrays

using Gridap
using Gridap.Helpers
using Gridap.CellValues: _test_iter_cell_value
using Gridap.CellValues: _test_index_cell_value

export CellArray
export CellMatrix
export CellVector

export IterCellArray
export IterCellMatrix
export IterCellVector

export IndexCellArray
export IndexCellMatrix
export IndexCellVector

export test_iter_cell_array
export test_index_cell_array

import Base: collect

# Iterable cell Arrays

const IterCellArray{T,N,A<:AbstractArray{T,N}} = IterCellValue{A}

const IterCellVector{T} = IterCellArray{T,1}

const IterCellMatrix{T} = IterCellArray{T,2}

# Indexable cell arrays

const IndexCellArray{T,N,A<:AbstractArray{T,N},D} = IndexCellValue{A,D}

const IndexCellVector{T} = IndexCellArray{T,1}

const IndexCellMatrix{T} = IndexCellArray{T,2}

# Cell Arrays

const CellArray{T,N} = Union{IterCellArray{T,N},IndexCellArray{T,N}}

const CellVector{T} = CellArray{T,1}

const CellMatrix{T} = CellArray{T,2}

collect(a::CellArray) = [ copy(ai) for ai in a ]

# Testers

function test_iter_cell_array(
  icv::CellArray{T,N},
  a::AbstractArray{<:AbstractArray{T,N}}) where {T,N}
  _test_iter_cell_value(icv,a)
end

function test_index_cell_array(
  icv::IndexCellArray{T,N},
  a::AbstractArray{<:AbstractArray{T,N}}) where {T,N}
  _test_index_cell_value(icv,a)
end

end # module

