
# Iterable cell fields

abstract type IterCellField{D,T} end

function iterate(::IterCellField{D,T})::Union{Nothing,Tuple{Field{D,T},Any}} where {D,T}
  @abstractmethod
end

function iterate(::IterCellField{D,T},state)::Union{Nothing,Tuple{Field{D,T},Any}} where {D,T}
  @abstractmethod
end

length(::IterCellField)::Int = @abstractmethod

eltype(::Type{C}) where C <: IterCellField{D,T} where {D,T} = Field{D,T}

IteratorEltype(::Type{C} where C <: IterCellField{D,T} where {D,T}) = EltypeUnknown()

# Indexable cell fields

abstract type IndexCellArray{T,N,A<:AbstractArray{T,N},D} <: AbstractArray{A,D} end

const CellField{D,T} = Union{IterCellField{D,T},IndexCellField{D,T}}

cellsize(self::CellArray,i::Int) = (s = cellsize(self); s[i])

celllength(self::CellArray) = prod(cellsize(self))
