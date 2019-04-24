module CellMaps

using Numa.Helpers

using Numa.Maps
using Numa.Maps: range_size
using Numa.Maps: MapFromUnaryOp
using Numa.Maps: MapFromBinaryOp
using Numa.Maps: FieldFromExpand
using Numa.Maps: FieldFromCompose
using Numa.Maps: FieldFromComposeExtended

using Numa.FieldValues
import Numa.FieldValues: inner, outer
using Numa.CellValues
using Numa.CellValues: CellArrayFromUnaryOp
using Numa.CellValues: CellArrayFromBroadcastUnaryOp
using Numa.CellValues: CellArrayFromBroadcastBinaryOp
using Numa.CellValues: CachedArray
using Numa.CellValues: setsize!

export CellMap
export CellField
export CellBasis
export CellGeomap

export CellMapValues
export CellBasisValues
export CellPoints

export ConstantCellMap

export expand
export varinner
export attachgeomap
export compose


import Base: +, -, *, /, ∘
import Base: iterate
import Base: length
import Base: eltype
import Base: size
import Base: getindex, setindex!

import Numa: evaluate, gradient
import Numa.CellValues: cellsize
import Numa.CellValues: inputcellarray, computesize, computevals!

# Iterable cell Maps
"""
Abstract object that traverses a set of cells and at every cell returns a
`Map{S,M,T,N}`
"""
abstract type IterCellMap{S,M,T,N} end
# @santiagobadia : Why don't put the result type R as template parameter,
# as for IndexCellMap ?

function iterate(::IterCellMap{S,M,T,N})::Union{Nothing,Tuple{Map{S,M,T,N},Any}} where {S,M,T,N}
  @abstractmethod
end

function iterate(::IterCellMap{S,M,T,N},state)::Union{Nothing,Tuple{Map{S,M,T,N},Any}} where {S,M,T,N}
  @abstractmethod
end

eltype(::Type{C}) where C <: IterCellMap{S,M,T,N} where {S,M,T,N} = Map{S,M,T,N}
# @santiagobadia :  I think this method must be overriden, better in CellMap?

# Indexable cell Maps

abstract type IndexCellMap{S,M,T,N,R<:Map{S,M,T,N}} <: AbstractVector{R} end

function getindex(::IndexCellMap{S,M,T,N,R}, ::Int)::R where {S,M,T,N,R}
  @abstractmethod
end

lastindex(x::IndexCellMap) = x[length(x)]

# Cell Maps
"""
Abstract object that for a given cell index returns a `Map{S,M,T,N}`
"""
const CellMap{S,M,T,N} = Union{IterCellMap{S,M,T,N},IndexCellMap{S,M,T,N}}
# santiagobadia : Problem if IterCellMap and IndexCellMap not same template types?
# Is this correct? IndexCellMap{S,M,T,N} when IndexCellMap{S,M,T,N,R}?

"""
Returns the evaluation of a `CellMap`
"""
function evaluate(::CellMap{S,M,T,N},::CellArray{S,M})::CellArray{T,N} where {S,M,T,N}
  @abstractmethod
end

"""
Returns another `CellMap` object that represents its gradient. Instances of `TG`
have a rank order a unit greater than the ones of `T`
"""
function gradient(::CellMap{S,M,T,N})::CellMap{S,M,TG,N} where {S,M,T<:FieldValue,N,TG}
  @abstractmethod
end

length(::CellMap)::Int = @abstractmethod

cellsize(::CellMap) = @abstractmethod
# @santiagobadia : What should I put here?

function Base.show(io::IO,self::CellMap)
  for (i, a) in enumerate(self)
    println(io,"$i -> $a")
  end
end

"""
Abstract type that represents a cell-wise basis for a field space,
where T is the type of value and D the dimension of the domain
"""
# const CellBasis{D,T} = CellMap{Point{D},1,T,2} where {D,T<:FieldValue}
const IterCellBasis{D,T} = IterCellMap{Point{D},1,T,2} where {D,T<:FieldValue}
const IndexCellBasis{D,T,R} = IndexCellMap{Point{D},1,T,2,R} where {D,T<:FieldValue,R}
const CellBasis{D,T} = Union{IterCellBasis{D,T},IndexCellBasis{D,T}}

"""
Abstract type that represents a cell-wise field, where
`T` stands for the type that represents the field at a point
(e.g., scalar, vector, tensor) and `D` stands for the space
dimension
"""
# const CellField{D,T} = CellMap{Point{D},1,T,1} where {D,T<:FieldValue}
const IterCellField{D,T} = IterCellMap{Point{D},1,T,1} where {D,T<:FieldValue}
const IndexCellField{D,T,R} = IndexCellMap{Point{D},1,T,1,R} where {D,T<:FieldValue,R}
const CellField{D,T} = Union{IterCellField{D,T},IndexCellField{D,T}}

"""
Abstract type representing a cell-wise transformation
between two geometrical domains
"""
const CellGeomap{D,Z} = CellField{D,Point{Z}}

# Abstract types for the input and output values
# of CellFields and CellBasis

"""
An array of points for each cell.
This type represent the objects where CellField and CellBasis are evaluated
"""
const CellPoints{D} = CellVector{Point{D}} where D

"""
Abstract type that represents a field with value of type T
evaluated at a collection of points in each cell
"""
const CellFieldValues{T} = CellVector{T} where T <: FieldValue

"""
Abstract type that represents a function basis with value of type T
evaluated at a collection of points in each cell
"""
const CellBasisValues{T} = CellArray{T,2} where T <: FieldValue

# @santiagobadia : Using the same as CellArray and CellValue

# Concrete structs

"""
Cell-wise map created from a `Map` of concrete type `R`
"""
struct ConstantCellMap{S,M,T,N,R <: Map{S,M,T,N}} <: IndexCellMap{S,M,T,N,R}
  map::R
  num_cells::Int
end

size(this::ConstantCellMap) = (this.num_cells,)

length(this::ConstantCellMap) = this.num_cells

"""
Evaluate a `ConstantCellMap` on a set of points represented with a
`CellArray{S,M}`
"""
function evaluate(self::ConstantCellMap{S,M,T,N,R},
  points::CellArray{S,M}) where {S,M,T,N,R}
  IterConstantCellMapValues(self.map,points)
end

"""
Computes the gradient of a `ConstantCellMap`
"""
function gradient(self::ConstantCellMap)
  gradfield = gradient(self.map)
  ConstantCellMap(gradfield,self.num_cells)
end

getindex(this::ConstantCellMap, i::Int) = this.map

firstindex(this::ConstantCellMap) = this.map

lastindex(this::ConstantCellMap) = this.map

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

function cellsize(this::IterConstantCellMapValues)
  return (range_size(this.map)..., cellsize(this.cellpoints)...)
end

@inline function Base.iterate(this::IterConstantCellMapValues{S,M,T,N,A,B}) where {S,M,T,N,A,B}
  # R = Base._return_type(evaluate,Tuple{A,B})
  # I could write a more general type by temp wrt R
  u = Array{T,N}(undef, cellsize(this))
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
  vsize = (range_size(this.map)..., size(a)...)
  # vsize = size(a)
  setsize!(v,vsize)
  evaluate!(this.map,a,v)
  state = (v, astate)
  (v, state)
end

const ConstantCellMapValues = IterConstantCellMapValues

include("Operators.jl")
include("CellBasisWithGeomap.jl")

end #module CellMaps
