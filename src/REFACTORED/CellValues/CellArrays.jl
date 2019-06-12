module CellArrays

using Gridap
using Gridap.Helpers

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

const IterCellArray = IterCellValue{A} where A<:AbstractArray

const IterCellVector = IterCellValue{A} where A<:AbstractVector

const IterCellMatrix = IterCellValue{A} where A<:AbstractMatrix

# Indexable cell arrays

const IndexCellArray = IndexCellValue{A,D} where {A<:AbstractArray,D}

const IndexCellVector = IndexCellValue{A,D} where {A<:AbstractVector,D}

const IndexCellMatrix = IndexCellValue{A,D} where {A<:AbstractMatrix,D}

# Cell Arrays

const CellArray = Union{IterCellArray{A},IndexCellArray{A,D}} where {A<:AbstractArray,D}

const CellVector = CellArray{A,D} where {A<:AbstractVector,D}

const CellMatrix = CellArray{A,D} where {A<:AbstractMatrix,D}

collect(a::CellArray) = [ copy(ai) for ai in a ]

# Testers

function test_iter_cell_array(
  icv::CellArray{<:AbstractArray{T,N}},
  a::AbstractArray{<:AbstractArray{T,N}}) where {T,N}
  test_iter_cell_value(icv,a)
end

function test_index_cell_array(
  icv::IndexCellArray{<:AbstractArray{T,N}},
  a::AbstractArray{<:AbstractArray{T,N}}) where {T,N}
  test_index_cell_value(icv,a)
end

end # module CellArrays

