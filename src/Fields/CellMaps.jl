module CellMaps

using Numa.Maps
using Numa.Maps: valsize

using Numa.Helpers
using Numa.CellValues
using Numa.CellValues: CachedArray
using Numa.CellValues: setsize!

import Base: iterate
import Base: length
import Base: eltype
import Base: size
import Base: getindex, setindex!

import Numa: evaluate, gradient
import Numa.CellValues: cellsize

# Iterable cell Maps

abstract type IterCellMap{S,M,T,N} end
# @santiagobadia : Why don't put the result type R as template parameter,
# as for IndexCellMap ?

function iterate(::IterCellMap{S,M,T,N})::Union{Nothing,Tuple{Map{S,M,T,N},Any}} where {S,M,T,N}
  @abstractmethod
end

function iterate(::IterCellMap{S,M,T,N},state)::Union{Nothing,Tuple{Map{S,M,T,N},Any}} where {S,M,T,N}
  @abstractmethod
end

length(::IterCellMap)::Int = @abstractmethod

eltype(::Type{C}) where C <: IterCellMap{S,M,T,N} where {S,M,T,N} = Map{S,M,T,N}

# Indexable cell Maps

abstract type IndexCellMap{S,M,T,N,R<:Map{S,M,T,N}} <: AbstractVector{R} end


function getindex(::IndexCellMap{S,M,T,N,R}, ::Int)::R where {S,M,T,N,R}
  @abstractmethod
end

lastindex(x::IndexCellMap) = x[length(x)]

# Cell Maps

const CellMap{S,M,T,N} = Union{IterCellMap{S,M,T,N},IndexCellMap{S,M,T,N}}
# santiagobadia : Problem if IterCellMap and IndexCellMap not same template types?

length(::CellMap)::Int = @abstractmethod

cellsize(::CellMap) = @abstractmethod
# @santiagobadia : What should I put here?



function Base.show(io::IO,self::CellMap)
  for (i, a) in enumerate(self)
    println(io,"$i -> $a")
  end
end
# @santiagobadia : Use the same as CellArray and CellValue

# Concrete structs

"""
Cell-wise map created from a `Map`
"""
struct ConstantCellMap{S,M,T,N,R} <: IndexCellMap{S,M,T,N,R}
  map::R
  num_cells::Int
end

size(this::ConstantCellMap) = (this.num_cells,)
length(this::ConstantCellMap) = this.num_cells

function ConstantCellMap(m::Map{S,M,T,N}, l::Int) where {S,M,T,N}
  R = typeof(m)
  ConstantCellMap{S,M,T,N,R}(m,l)
end

function evaluate(self::ConstantCellMap{S,M,T,N,R},points::CellArray{S,M}) where {S,M,T,N,R,P}
  IterConstantCellMapValues(self.map,points)
end

function gradient(self::ConstantCellMap)
  gradfield = gradient(self.field)
  ConstantCellMap(gradfield)
end

getindex(this::ConstantCellMap, i::Int) = this.map
firstindex(this::ConstantCellMap) = this.map
lastindex(this::ConstantCellMap) = this.map

# CellMapValues

struct IterConstantCellMapValues{S,M,T,N,A<:Map{S,M,T,N},B<:CellArray{S,M}} <: IterCellArray{T,N}
  map::A
  cellpoints::B
end

cellsize(this::IterConstantCellMapValues) = cellsize(this.cellpoints)

@inline function Base.iterate(this::IterConstantCellMapValues{S,M,T,N,A,B}) where {S,M,T,N,A,B}
  R = Base._return_type(evaluate,Tuple{A,B})
  # @santiagobadia : Here is the problem...
  #  a field should be S,0,T,0 and after evaluation, it would take e.g., S,1
  # and return T,1... i.e. T,N+1
  u = Array{T,N}(undef,cellsize(this.cellpoints))
  v = CachedArray(u)
  anext = iterate(this.cellpoints)
  if anext === nothing; return nothing end
  iteratekernel(this,anext,v)
end

@inline function Base.iterate(this::IterConstantCellMapValues{S,M,T,N,A,B},state) where {S,M,T,N,A,B}
  v, astate = state
  anext = iterate(this.cellpoints,astate)
  if anext === nothing; return nothing end
  iteratekernel(this,anext,v)
end

function iteratekernel(this::IterConstantCellMapValues,next,v)
  a, astate = next
  # vsize = computesize(size(a))
  vsize = size(a)
  setsize!(v,vsize)
  evaluate!(this.map,a,v)
  # computevals!(this,a,v)
  state = (v, astate)
  (v, state)
  # @santiagobadia : I don't understand the last step, I have copied from a
  # similar situation in other part of the code, but
end

const ConstantCellMapValues = IterConstantCellMapValues

end #module CellMaps
