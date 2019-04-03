
abstract type IterCellValue{T} end

function iterate(::IterCellValue{T})::Union{Nothing,Tuple{T,Any}} where T
  @abstractmethod
end

function iterate(::IterCellValue{T},state)::Union{Nothing,Tuple{T,Any}} where T
  @abstractmethod
end

length(::IterCellValue)::Int = @abstractmethod

eltype(::Type{C}) where C <: IterCellValue{T} where T = T

function Base.show(io::IO,self::IterCellValue)
  for (i, a) in enumerate(self)
    println(io,"$i -> $a")
  end
end

abstract type IndexCellValue{T} <: AbstractVector{T} end

const CellValue{T} = Union{IterCellValue{T},IndexCellValue{T}}

const CellArray{T,N} = CellValue{<:AbstractArray{T,N}} where {T,N}

const CellVector{T} = CellArray{T,1} where T

function cellsize(::CellArray{T,N})::NTuple{N,Int}  where {T,N}
  @abstractmethod
end

cellsize(self::CellArray,i::Int) = (s = cellsize(self); s[i])

celllength(self::CellArray) = prod(cellsize(self))

