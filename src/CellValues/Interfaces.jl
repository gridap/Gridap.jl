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

# @santiagobadia : Check all the methods to be defined for abstract iterators
# and arrays

# Indexable cell Values

abstract type IndexCellValue{T,N} <: AbstractArray{T,N} end
# @santiagobadia : I think it should be a sub-type of AbstractArray{T,N}
# since the arrays can have variable size at different cells

function getindex(::IndexCellValue{T,N}, ::Int)::R where {T,N,R}
  @abstractmethod
end

# function getindex(::IndexCellValue{T,N}, ::Vararg{Int,D})::R where {T,N,R,D}
  # @abstractmethod
# end
# @santiagobadia : If I uncomment this method concrete implementations with
# IndexLinear do not work, because they get here (?)

size(x::IndexCellValue) = @abstractmethod

lastindex(x::IndexCellValue) = x[length(x)]
# @santiagobadia : Not true if it is sub-typing AbsttractArray{T,N}

IndexStyle(::IndexCellValue) = IndexLinear()

# Cell Values

const CellValue{T} = Union{IterCellValue{T},IndexCellValue{T}}

cellsize(::CellValue) = ()
# @santiagobadia : Does it have sense here?

# Iterable cell Arrays

const IterCellArray{T,N} = IterCellValue{AbstractArray{T,N}}

IteratorEltype(::Type{C} where C <: IterCellArray{T,N} where {T,N}) = EltypeUnknown()

# Indexable cell arrays

const IndexCellArray{T,N,A<:AbstractArray{T,N},D} = IndexCellValue{A,D}

# const IndexCellVector{T,A,D} = IndexCellArray{T,1,A,D}
# @santiagobadia : I don't get why
# IndexCellVector{T} <: IndexCellArray{T,1} not true????

# Cell Arrays

const CellArray{T,N} = Union{IterCellArray{T,N},IndexCellArray{T,N}}
# @santiagobadia : Related to the previous comment, it is not valid to sub-type
# CellArray{T,N}. Same for CellField and CellBasis

const CellVector{T} = CellArray{T,1} where T

const CellMatrix{T} = CellArray{T,2} where T

cellsize(self::CellArray,i::Int) = (s = cellsize(self); s[i])

celllength(self::CellArray) = prod(cellsize(self))

function Base.show(io::IO,self::CellValue)
  for (i, a) in enumerate(self)
    println(io,"$i -> $a")
  end
end

const IterData{T} = Union{CellValue{T},AbstractArray{T}}
