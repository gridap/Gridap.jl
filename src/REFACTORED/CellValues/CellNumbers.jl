module CellNumbers

using Gridap
using Gridap.Helpers

using StaticArrays

export NumberLike
export CellNumber
export IterCellNumber
export IndexCellNumber

export test_index_cell_number
export test_iter_cell_number

const NumberLike = Union{Number,SArray}

# Iterable cell Numbers

const IterCellNumber = IterCellValue{A} where A<:NumberLike

# Indexable cell arrays

const IndexCellNumber = IndexCellValue{A,D} where {A<:NumberLike,D}

# Cell Numbers

const CellNumber = Union{IterCellNumber{A},IndexCellNumber{A,D}} where {A<:NumberLike,D}

# Testers

function test_iter_cell_number(icv::CellNumber, a::AbstractArray{<:NumberLike})
  test_iter_cell_value(icv,a)
end

function test_index_cell_number(icv::IndexCellNumber, a::AbstractArray{<:NumberLike})
  test_index_cell_value(icv,a)
end

end # module CellNumbers
