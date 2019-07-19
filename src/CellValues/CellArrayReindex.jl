module CellArrayReindex

using Gridap
using Gridap.CellValuesGallery
using Gridap.CachedArrays
using Base: @pure

import Gridap: reindex
import Base: iterate
import Base: length
import Base: size
import Base: getindex
import Base: IndexStyle

function reindex(values::IndexCellArray, indices::CellValue{<:IndexLike})
  IterCellArrayFromReindex(values,indices)
end

function reindex(values::IndexCellArray, indices::IndexCellValue{<:IndexLike})
  IndexCellArrayFromReindex(values,indices)
end

struct IterCellArrayFromReindex{T,N,V,I} <: IterCellValue{CachedArray{T,N,Array{T,N}}}
  values::V
  indices::I
end

function IterCellArrayFromReindex(
  values::IndexCellArray{T,N}, indices::CellValue{<:IndexLike}) where {T,N}
  V = typeof(values)
  I = typeof(indices)
  IterCellArrayFromReindex{T,N,V,I}(values,indices)
end

@inline function iterate(s::IterCellArrayFromReindex{T,N}) where {T,N}
  inext = iterate(s.indices)
  cache = CachedArray(T,N)
  _iterate(s,inext,cache)
end

@inline function iterate(s::IterCellArrayFromReindex,state)
  istate, cache = state
  inext = iterate(s.indices,istate)
  _iterate(s,inext,cache)
end

@inline function _iterate(s,inext,cache)
  if inext === nothing; return nothing end
  i, istate = inext
  v = s.values[i]
  _update_cache!(cache,v)
  state = (istate, cache)
  (cache, state)
end

function _update_cache!(cache,v)
  setsize!(cache,size(v))
  for (j,vi) in zip(eachindex(cache),v)
    cache[j] = vi
  end
end

length(s::IterCellArrayFromReindex) = length(s.indices)

struct IndexCellArrayFromReindex{T,N,D,V,I} <: IndexCellValue{CachedArray{T,N,Array{T,N}},D}
  values::V
  indices::I
  cache::CachedArray{T,N,Array{T,N}}
end

function IndexCellArrayFromReindex(
  values::IndexCellArray{T,N}, indices::IndexCellValue{<:IndexLike,D}) where {D,T,N}
  V = typeof(values)
  I = typeof(indices)
  cache = CachedArray(T,N)
  IndexCellArrayFromReindex{T,N,D,V,I}(values,indices,cache)
end

function getindex(
  s::IndexCellArrayFromReindex{T,D}, i::Vararg{<:Integer,D}) where {T,D}
  j = s.indices[i...]
  v = s.values[j]
  _update_cache!(s.cache,v)
  s.cache
end

size(s::IndexCellArrayFromReindex) = size(s.indices)

@pure function IndexStyle(
  ::Type{<:IndexCellArrayFromReindex{T,N,D,V,I}}) where {T,N,D,V,I}
  IndexStyle(I)
end

end # module
