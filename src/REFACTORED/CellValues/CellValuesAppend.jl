module CellValuesAppend

using Base: @propagate_inbounds

using Gridap
using Gridap.Helpers
using Gridap.CachedArrays

export append
export local_append

import Base: iterate
import Base: length
import Base: eltype
import Base: size
import Base: getindex
import Base: IndexStyle

function append(cvs::Vararg{IndexCellValue{T,1}}) where T
  IndexCellValueByGlobalAppend(cvs...)
end

function local_append(offset::Int,a::IndexCellVector,b::IndexCellVector)
  IndexCellVectorByLocalAppend(offset,a,b)
end

function local_append(offset::NTuple{M,Int},a::Vararg{<:IndexCellVector}) where M
  IndexCellVectorByLocalAppend(offset,a...)
end

struct IndexCellValueByGlobalAppend{
  T,V<:IndexCellValue{T,1},W<:IndexCellValue{T,1}} <: IndexCellValue{T,1}
  cvs1::V
  cvs2::W
  offset::Int
end

function IndexCellValueByGlobalAppend(cvs1::IndexCellValue{T,1}, cvs2::IndexCellValue{T,1}) where {T}
  offset = length(cvs1)
  IndexCellValueByGlobalAppend(cvs1, cvs2, offset)
end

function IndexCellValueByGlobalAppend(cvs::Vararg{IndexCellValue{T,1},N}) where {T,N}
  aux = cvs[1]
  for i in 2:N
    aux = IndexCellValueByGlobalAppend(aux, cvs[i])
  end
  return aux
end

@propagate_inbounds function getindex(self::IndexCellValueByGlobalAppend,cell::Int)
  if cell <= self.offset
    return self.cvs1[cell]
  else
    return self.cvs2[cell-self.offset]
  end
end

size(self::IndexCellValueByGlobalAppend) = (length(self.cvs1) + length(self.cvs2),)

struct IndexCellVectorByLocalAppend{
  T,V<:IndexCellVector{T},W<:IndexCellVector{T}} <: IndexCellVector{T,CachedArray{T,1,Array{T,1}},1}
  offset::Int
  cvs1::V
  cvs2::W
  cache::CachedArray{T,1,Array{T,1}}
end

function IndexCellVectorByLocalAppend(
  offset::Int,cvs1::IndexCellVector{T}, cvs2::IndexCellVector{T}) where T
  @assert length(cvs1) == length(cvs2)
  cache = CachedArray(T,1)
  IndexCellVectorByLocalAppend(offset, cvs1, cvs2, cache)
end

function IndexCellVectorByLocalAppend(
  offset_vec::NTuple{M,Int}, cvs::Vararg{<:IndexCellVector{T},N}) where {T,N,M}
  @assert M == N-1
  aux = cvs[1]
  for i in 2:N
    aux = IndexCellVectorByLocalAppend(offset_vec[i-1], aux, cvs[i])
  end
  return aux
end

@propagate_inbounds function getindex(self::IndexCellVectorByLocalAppend,cell::Integer)
  v1 = self.cvs1[cell]
  v2 = self.cvs2[cell]
  l1 = length(v1)
  l2 = length(v2)
  l = l1 + l2
  s = (l,)
  setsize!(self.cache,s)
  for i in 1:l1
    @inbounds self.cache[i] = v1[i]
  end
  for i in 1:l2
    @inbounds self.cache[i+l1] = v2[i] + self.offset
  end
  self.cache
end

length(self::IndexCellVectorByLocalAppend) = length(self.cvs1)

size(self::IndexCellVectorByLocalAppend) = (length(self),)

end # module
