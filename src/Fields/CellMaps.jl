module CellMaps

using Numa.Maps
using Numa.Helpers

import Base: iterate
import Base: length
import Base: eltype
import Base: size
import Base: getindex, setindex!

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

cellsize(::CellMap) = ()
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
  length::Int
end

function ConstantCellMap(m::Map{S,M,T,N}, l::Int) where {S,M,T,N}
  R = typeof(m)
  ConstantCellMap{S,M,T,N,R}(m,l)
end

function evaluate(self::ConstantCellMap{S,M,T,N,R},points::AbstractVector{P}) where {S,M,T,N,R,P}
  ConstantCellMapValues(self.maps,points)
end
# @santiagobadia : Is the P type OK? I would put <: AbstractArray{S,M} but I guess
# I will have problems for M=0, i.e., just a value of S...

function gradient(self::ConstantCellMap)
  gradfield = gradient(self.field)
  ConstantCellMap(gradfield)
end

getindex(this::ConstantCellMap, i::Int) = this.map
firstindex(this::ConstantCellMap) = this.map
lastindex(this::ConstantCellMap) = this.map

end #module CellMaps
