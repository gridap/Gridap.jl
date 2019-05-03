"""
Cell-wise map created from a `Map` of concrete type `R`
"""
const ConstantCellMap{S,M,T,N,R<:Map{S,M,T,N}} = ConstantCellValue{R}

function ConstantCellMap(a::Map,l::Int)
  return ConstantCellValue(a,l)
end
"""
Evaluate a `ConstantCellMap` on a set of points represented with a
`CellArray{S,M}`
"""
function evaluate(self::ConstantCellMap{S,M,T,N,R},
  points::CellArray{S,M}) where {S,M,T,N,R}
  IterConstantCellMapValues(self.value,points)
end

"""
Computes the gradient of a `ConstantCellMap`
"""
function gradient(self::ConstantCellMap)
  gradfield = gradient(self.value)
  ConstantCellMap(gradfield,self.length)
end

getindex(this::ConstantCellMap, i::Int) = this.value

# firstindex(this::ConstantCellMap) = this.value

# lastindex(this::ConstantCellMap) = this.value

# CellMapValues

"""
Object that represents the (lazy) evaluation of a `CellMap{S,M,T,N}` on a set
of points represented with a `CellArray{S,M,T,N}`. The result is a sub-type of
`IterCellArray{T,N}`. Its template parameters are `{S,M,T,N,A,B}`, where `A`
stands for the concrete sub-type of `Map{S,M,T,N}` and `B` stands for the
concrete sub-type of `CellArray{S,M}`
"""
struct IterConstantCellMapValues{S,M,T,N,A<:Map{S,M,T,N},B<:CellArray{S,M}} <: IterCellArray{T,N}
  map::A
  cellpoints::B
end

length(this::IterConstantCellMapValues) = length(this.map)

function cellsize(this::IterConstantCellMapValues)
  return_size(this.map, cellsize(this.cellpoints))
  # @notimplemented
end

@inline function Base.iterate(this::IterConstantCellMapValues{S,M,T,N,A,B}) where {S,M,T,N,A,B}
  u = Array{T,N}(undef, cellsize(this)...)
  v = CachedArray(u)
  anext = iterate(this.cellpoints)
  if anext === nothing; return nothing end
  iteratekernel(this,anext,v)
end
# @santiagobadia :  Create version with ConstantCellArray{Point{D},1}
# Does it pay the price? Performance?

@inline function Base.iterate(this::IterConstantCellMapValues{S,M,T,N,A,B},state) where {S,M,T,N,A,B}
  v, astate = state
  anext = iterate(this.cellpoints,astate)
  if anext === nothing; return nothing end
  iteratekernel(this,anext,v)
end

function iteratekernel(this::IterConstantCellMapValues,next,v)
  a, astate = next
  vsize = return_size(this.map,size(a))
  setsize!(v,vsize)
  evaluate!(this.map,a,v)
  state = (v, astate)
  (v, state)
end

const ConstantCellMapValues = IterConstantCellMapValues
