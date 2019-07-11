module CellValuesReindex

using Gridap
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
  IterCellValueFromReindex(values,indices)
end

function reindex(values::IndexCellValue, indices::IndexCellValue{<:IndexLike})
  IndexCellValueFromReindex(values,indices)
end

struct IterCellValueFromReindex{T,V,I} <: IterCellValue{T}
  values::V
  indices::I
end

function IterCellValueFromReindex(
  values::IndexCellValue{T}, indices::CellValue{<:IndexLike}) where T
  V = typeof(values)
  I = typeof(indices)
  IterCellValueFromReindex{T,V,I}(values,indices)
end

function iterate(s::IterCellValueFromReindex)
  inext = iterate(s.indices)
  _iterate(s,inext)
end

function iterate(s::IterCellValueFromReindex,state)
  inext = iterate(s.indices,state)
  _iterate(s,inext)
end

function _iterate(s,inext)
  if inext === nothing; return nothing end
  i, istate = inext
  (s.values[i], istate)
end

length(s::IterCellValueFromReindex) = length(s.indices)

struct IndexCellValueFromReindex{T,D,V,I} <: IndexCellValue{T,D}
  values::V
  indices::I
end

function IndexCellValueFromReindex(
  values::IndexCellValue{T}, indices::IndexCellValue{<:IndexLike,D}) where {D,T}
  V = typeof(values)
  I = typeof(indices)
  IndexCellValueFromReindex{T,D,V,I}(values,indices)
end

function getindex(
  s::IndexCellValueFromReindex{T,D}, i::Vararg{<:Integer,D}) where {T,D}
  j = s.indices[i...]
  s.values[j]
end

size(s::IndexCellValueFromReindex) = size(s.indices)

@pure function IndexStyle(
  ::Type{<:IndexCellValueFromReindex{T,D,V,I}}) where {T,D,V,I}
  IndexStyle(I)
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
