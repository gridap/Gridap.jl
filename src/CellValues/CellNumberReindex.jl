module CellNumberReindex

using Gridap
using Gridap.CellValuesGallery
using Gridap.CellMapApply: CellMapFromKernel
using Gridap.CellMapApply: IndexCellMapFromKernel
using Base: @pure

import Gridap: reindex
import Base: iterate
import Base: length
import Base: size
import Base: getindex
import Base: IndexStyle

function reindex(values::IndexCellNumber, indices::CellValue{<:IndexLike})
  IterCellNumberFromReindex(values,indices)
end

function reindex(values::IndexCellNumber, indices::IndexCellValue{<:IndexLike})
  IndexCellNumberFromReindex(values,indices)
end

struct IterCellNumberFromReindex{T,V,I} <: IterCellNumber{T}
  values::V
  indices::I
end

function IterCellNumberFromReindex(
  values::IndexCellNumber{T}, indices::CellValue{<:IndexLike}) where T
  V = typeof(values)
  I = typeof(indices)
  IterCellNumberFromReindex{T,V,I}(values,indices)
end

function iterate(s::IterCellNumberFromReindex)
  inext = iterate(s.indices)
  _iterate(s,inext)
end

function iterate(s::IterCellNumberFromReindex,state)
  inext = iterate(s.indices,state)
  _iterate(s,inext)
end

function _iterate(s,inext)
  if inext === nothing; return nothing end
  i, istate = inext
  (s.values[i], istate)
end

length(s::IterCellNumberFromReindex) = length(s.indices)

struct IndexCellNumberFromReindex{T,D,V,I} <: IndexCellNumber{T,D}
  values::V
  indices::I
end

function IndexCellNumberFromReindex(
  values::IndexCellNumber{T}, indices::IndexCellValue{<:IndexLike,D}) where {D,T}
  V = typeof(values)
  I = typeof(indices)
  IndexCellNumberFromReindex{T,D,V,I}(values,indices)
end

function getindex(
  s::IndexCellNumberFromReindex{T,D}, i::Vararg{<:Integer,D}) where {T,D}
  j = s.indices[i...]
  s.values[j]
end

size(s::IndexCellNumberFromReindex) = size(s.indices)

@pure function IndexStyle(
  ::Type{<:IndexCellNumberFromReindex{T,D,V,I}}) where {T,D,V,I}
  IndexStyle(I)
end

end # module
