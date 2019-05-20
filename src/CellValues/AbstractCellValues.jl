module AbstractCellValues

using Gridap.Helpers
using StaticArrays

export CellValue
export CellArray
export CellMatrix
export CellVector

export IterCellValue
export IterCellArray
export IterCellMatrix
export IterCellVector

export IndexCellValue
export IndexCellArray
export IndexCellMatrix
export IndexCellVector

import Base: iterate
import Base: length
import Base: eltype
import Base: size
import Base: getindex
import Base: IndexStyle
import Base: show
import Base: collect

export cellsize
export celllength

# Iterable cell Values

abstract type IterCellValue{T} end

function iterate(::IterCellValue{T})::Union{Nothing,Tuple{T,Any}} where T
  @abstractmethod
end

function iterate(::IterCellValue{T},state)::Union{Nothing,Tuple{T,Any}} where T
  @abstractmethod
end

length(::IterCellValue)::Int = @abstractmethod

eltype(::Type{C}) where C <: IterCellValue{T} where T = T

# Indexable cell Values

abstract type IndexCellValue{T,N} <: AbstractArray{T,N} end

function getindex(::IndexCellValue{T,N}, ::Vararg{Int,N})::T where {T,N}
  @abstractmethod
end

size(x::IndexCellValue) = @abstractmethod

IndexStyle(::Type{<:IndexCellValue{T,N}} where {T,N}) = IndexLinear()

# Cell Values

const CellValue{T} = Union{IterCellValue{T},IndexCellValue{T}}

# Iterable cell Arrays

const IterCellArray{T,N,A<:AbstractArray{T,N}} = IterCellValue{A}

const IterCellVector{T,A} = IterCellArray{T,1,A}

const IterCellMatrix{T,A} = IterCellArray{T,2,A}

function cellsize(::IterCellArray{T,N})::NTuple{N,Int} where {T,N}
  @abstractmethod
end

# Indexable cell arrays

const IndexCellArray{T,N,A<:AbstractArray{T,N},D} = IndexCellValue{A,D}

const IndexCellVector{T,A<:AbstractArray{T,1},D} = IndexCellArray{T,1,A,D}

const IndexCellMatrix{T,A<:AbstractArray{T,2},D} = IndexCellArray{T,2,A,D}

function cellsize(::IndexCellArray{T,N})::NTuple{N,Int} where {T,N}
  @abstractmethod
end

# Cell Arrays

const CellArray{T,N} = Union{IterCellArray{T,N},IndexCellArray{T,N}}

const CellVector{T} = CellArray{T,1}

const CellMatrix{T} = CellArray{T,2}


# Misc. methods that can be implemented at this level

cellsize(self::CellValue) = ()

cellsize(self::CellValue{<:SArray}) = ()

cellsize(self::CellArray,i::Int) = (s = cellsize(self); s[i])

celllength(self::CellArray) = prod(cellsize(self))

function show(io::IO,self::CellValue)
  for (i, a) in enumerate(self)
    println(io,"$i -> $a")
  end
end

collect(a::CellArray) = [ copy(ai) for ai in a ]

end # module AbstractCellValues
