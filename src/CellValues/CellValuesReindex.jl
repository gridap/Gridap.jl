module CellValuesReindex

using Gridap
using Gridap.Helpers
using Gridap.CellValuesGallery
using Gridap.CellMapApply: CellMapFromKernel
using Gridap.CellMapApply: IndexCellMapFromKernel
using Base: @pure

export reindex
export IndexLike

import Base: iterate
import Base: length
import Base: size
import Base: getindex
import Base: IndexStyle

const IndexLike = Union{Integer,CartesianIndex}

function reindex(values::IndexCellValue, indices::Vector{<:IndexLike})
  i = CellValueFromArray(indices)
  reindex(values,i)
end

function reindex(values::IndexCellValue, indices::CellValue{<:IndexLike})
  @notimplemented
end

function reindex(values::IndexCellValue, indices::IndexCellValue{<:IndexLike})
  @notimplemented
end

# Efficient implementations of reindex for concrete types

function reindex(cm::CellMapFromKernel,indices::CellValue{<:IndexLike})
  _reindex_cmfk(cm,indices)
end

function reindex(cm::IndexCellMapFromKernel,indices::CellValue{<:IndexLike})
  _reindex_cmfk(cm,indices)
end

function reindex(cm::IndexCellMapFromKernel,indices::IndexCellValue{<:IndexLike})
  _reindex_cmfk(cm,indices)
end

function _reindex_cmfk(cm,indices)
  k = cm.kernel
  rs = [ reindex(mi,indices) for mi in cm.cellvalues ]
  apply(k,rs...)
end

end # module
