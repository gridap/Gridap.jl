module Append

using Base: @propagate_inbounds
using StaticArrays

using Gridap.Helpers
using Gridap.CellValues
using Gridap.CachedArrays

export IndexCellValueByGlobalAppend
export IndexCellValueByLocalAppend
export IndexCellValueByLocalAppendWithOffset

import Base: iterate
import Base: length
import Base: eltype
import Base: size
import Base: getindex
import Base: IndexStyle
import Gridap.CellValues: cellsize

#TODO test the structs in this module with the tester function `test_index_cell_value` and `test_index_cell_array`
# depending of the case

struct IndexCellValueByGlobalAppend{T,V<:IndexCellValue{T,1},W<:IndexCellValue{T,1}} <: IndexCellValue{T,1}
  cvs1::V
  cvs2::W
  offset::Int
end

@propagate_inbounds function getindex(self::IndexCellValueByGlobalAppend,cell::Int)
  if cell <= self.offset
    return self.cvs1[cell]
  else
    return self.cvs2[cell-self.offset]
  end
end

size(self::IndexCellValueByGlobalAppend) = (length(self.cvs1) + length(self.cvs2),)

IndexStyle(::Type{IndexCellValueByGlobalAppend{T,N,V}}) where {T,N,V} = IndexStyle(V)

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

function cellsize(self::IndexCellValueByGlobalAppend)
  max( cellsize(self.cvs1), cellsize(self.cvs2))
end

# TODO I think we need to inherit from IndexCellVector (see below)
struct IndexCellValueByLocalAppend{T,V<:IndexCellValue{T,1},W<:IndexCellValue{T,1}} <: IndexCellValue{T,1}
  offset::Int
  cvs1::V
  cvs2::W
end

function IndexCellValueByLocalAppend(cvs1::IndexCellValue{T,1}, cvs2::IndexCellValue{T,1}) where {T}
  IndexCellValueByLocalAppend(0, cvs1, cvs2)
end

function IndexCellValueByLocalAppendWithOffset(cvs1::IndexCellValue{T,1}, cvs2::IndexCellValue{T,1}) where {T}
  m = max(cvs1[end]...)
  for i in cvs1
    m = max(i...,m)
  end
  IndexCellValueByLocalAppend(m, cvs1, cvs2)
end

@propagate_inbounds function getindex(self::IndexCellValueByLocalAppend,cell::Int)
  return [ self.cvs1[cell]..., (self.cvs2[cell].+self.offset)...] #TODO this creates a new object. Use CachedVector instead
  #TODO We are returning a vector, thus we need to inherit from IndexCellVector
end

size(self::IndexCellValueByLocalAppend) = size(self.cvs1)

function cellsize(self::IndexCellValueByLocalAppend)
  cs1 = cellsize(self.cvs1)
  cs2 = cellsize(self.cvs2)
  tuple(  (s1+s2 for (s1,s2) in zip(cs1,cs2))... )
end

# IndexStyle(::Type{IndexCellValueByLocalAppend{T,N,V}}) where {T,N,V} =  IndexLinear()

function IndexCellValueByLocalAppend(cvs::Vararg{IndexCellValue{T,1},N}) where {T,N}
  aux = cvs[1]
  for i in 2:N
    aux = IndexCellValueByLocalAppend(0, aux, cvs[i])
  end
  return aux
end

function IndexCellValueByLocalAppendWithOffset(offset_vec::NTuple{M,Int}, cvs::Vararg{IndexCellValue{T,1},N}) where {T,N,M}
  @assert M >= N-1
  aux = cvs[1]
  for i in 2:N
    aux = IndexCellValueByLocalAppend(offset_vec[i-1], aux, cvs[i])
  end
  return aux
end

function IndexCellValueByLocalAppendWithOffset(cvs::Vararg{IndexCellValue{T,1},N}) where {T,N}
  aux = cvs[1]
  for i in 2:N
    aux = IndexCellValueByLocalAppendWithOffset(aux, cvs[i])
  end
  return aux
end

end # module Append
