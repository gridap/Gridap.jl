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

abstract type IndexCellValue{T} <: AbstractVector{T} end

# Cell Values

const CellValue{T} = Union{IterCellValue{T},IndexCellValue{T}}

# Iterable cell Arrays

abstract type CellArray{T,N} end

const CellVector{T} = CellArray{T,1} where T

const CellMatrix{T} = CellArray{T,2} where T

function iterate(::CellArray{T,N})::Union{Nothing,Tuple{AbstractArray{T,N},Any}} where {T,N}
  @abstractmethod
end

function iterate(::CellArray{T,N},state)::Union{Nothing,Tuple{AbstractArray{T,N},Any}} where {T,N}
  @abstractmethod
end

length(::CellArray)::Int = @abstractmethod

IteratorEltype(::Type{C} where C <: CellArray{T,N} where {T,N}) = EltypeUnknown()

cellsize(self::CellArray,i::Int) = (s = cellsize(self); s[i])

celllength(self::CellArray) = prod(cellsize(self))

# Indexable cell Arrays
# We don't extend from AbstractVector, since in general we do not know
# the type of the returned matrix

abstract type IndexCellArray{T,N} <: CellArray{T,N} end

function getindex(::IndexCellArray{T,N},cell::Int)::AbstractArray{T,N} where {T,N}
  @abstractmethod
end

@inline iterate(self::IndexCellArray) = iterate(self,0)

@inline function iterate(self::IndexCellArray,state::Int)
  if length(self) == state
    nothing
  else
    k = state+1
    (self[k],k)
  end
end

# Cell Data

const CellData{T} = Union{CellValue{T},CellArray{T}}

function Base.show(io::IO,self::CellData)
  for (i, a) in enumerate(self)
    println(io,"$i -> $a")
  end
end

