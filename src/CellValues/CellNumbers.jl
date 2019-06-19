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

const IterCellNumber{A<:NumberLike} = IterCellValue{A}

# Indexable cell arrays

const IndexCellNumber{A<:NumberLike,D} = IndexCellValue{A,D}

# Cell Numbers

const CellNumber{A<:NumberLike} = Union{IterCellNumber{A},IndexCellNumber{A}}

# Testers

function test_iter_cell_number(icv::CellNumber, a::AbstractArray{<:NumberLike})
  test_iter_cell_value(icv,a)
end

function test_index_cell_number(icv::IndexCellNumber, a::AbstractArray{<:NumberLike})
  test_index_cell_value(icv,a)
end

end # module CellNumbers
