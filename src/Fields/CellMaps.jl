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

const CellMap{S,M,T,N,R} = Union{IterCellMap{S,M,T,N},IndexCellMap{S,M,T,N,R}}
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

end #module CellMaps
