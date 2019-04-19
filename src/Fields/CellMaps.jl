module CellMaps

using Numa.Helpers

# Iterable cell Maps

abstract type IterCellMap{T} end

function iterate(::IterCellMap{T})::Union{Nothing,Tuple{T,Any}} where T
  @abstractmethod
end

function iterate(::IterCellMap{T},state)::Union{Nothing,Tuple{T,Any}} where T
  @abstractmethod
end

length(::IterCellMap)::Int = @abstractmethod

eltype(::Type{C}) where C <: IterCellMap{T} where T = T

# Indexable cell Maps

abstract type IndexCellMap{T,N} <: AbstractArray{T,N} end

# Cell Maps

const CellMap{T} = Union{IterCellMap{T},IndexCellMap{T}}

cellsize(::CellMap) = ()
# @santiagobadia : What should I put here?

function Base.show(io::IO,self::CellMap)
  for (i, a) in enumerate(self)
    println(io,"$i -> $a")
  end
end
# @santiagobadia : Use the same as CellArray and CellValue

end #module CellMaps
