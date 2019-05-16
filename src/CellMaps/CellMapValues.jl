module CellMapValues

using Gridap
using Gridap.CellMaps
using Gridap.CellValues
using Gridap.CachedArrays

export CellMapValue

import Base: iterate, length
import Gridap.CellValues: cellsize

struct CellMapValue{
  S,M,T,N,A<:CellMap{S,M,T,N},B<:CellArray{S,M}} <: IterCellArray{T,N,CachedArray{T,N,Array{T,N}}}
  cellmap::A
  cellarray::B
end

@inline function iterate(self::CellMapValue{S,M,T,N}) where {S,M,T,N}
  u = Array{T,N}(undef, cellsize(self))
  v = CachedArray(u)
  anext = iterate(self.cellmap)
  bnext = iterate(self.cellarray)
  if anext === nothing; return nothing end
  if bnext === nothing; return nothing end
  iteratekernel(self,anext,bnext,v)
end

@inline function iterate(self::CellMapValue,state)
  v, astate, bstate = state
  anext = iterate(self.cellmap,astate)
  bnext = iterate(self.cellarray,bstate)
  if anext === nothing; return nothing end
  if bnext === nothing; return nothing end
  iteratekernel(self,anext,bnext,v)
end

function iteratekernel(self::CellMapValue,anext,bnext,v)
  a, astate = anext
  b, bstate = bnext
  vsize = return_size(a,size(b))
  setsize!(v,vsize)
  evaluate!(a,b,v)
  state = (v, astate, bstate)
  (v,state)
end

function length(self::CellMapValue)
  @assert length(self.cellmap) == length(self.cellarray)
  length(self.cellmap)
end

function cellsize(self::CellMapValue) where {T,N}
  return_size(self.cellmap,cellsize(self.cellarray))
end

end # module CellMapValues
