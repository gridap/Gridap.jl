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

# Cell Values

const CellValue{T} = Union{IterCellValue{T},IndexCellValue{T}}

cellsize(::CellValue) = ()

# Iterable cell Arrays

abstract type IterCellArray{T,N} end

function iterate(::IterCellArray{T,N})::Union{Nothing,Tuple{AbstractArray{T,N},Any}} where {T,N}
  @abstractmethod
end

function iterate(::IterCellArray{T,N},state)::Union{Nothing,Tuple{AbstractArray{T,N},Any}} where {T,N}
  @abstractmethod
end

length(::IterCellArray)::Int = @abstractmethod

IteratorEltype(::Type{C} where C <: IterCellArray{T,N} where {T,N}) = EltypeUnknown()

# Indexable cell arrays

abstract type IndexCellArray{T,N,A<:AbstractArray{T,N},D} <: AbstractArray{A,D} end

# Cell Arrays

const CellArray{T,N} = Union{IterCellArray{T,N},IndexCellArray{T,N}}

const CellVector{T} = CellArray{T,1} where T

const CellMatrix{T} = CellArray{T,2} where T

cellsize(self::CellArray,i::Int) = (s = cellsize(self); s[i])

celllength(self::CellArray) = prod(cellsize(self))

# Cell Data

const CellData{T} = Union{CellValue{T},CellArray{T}}

function Base.show(io::IO,self::CellData)
  for (i, a) in enumerate(self)
    println(io,"$i -> $a")
  end
end

cellsize(::CellData) = @abstractmethod
