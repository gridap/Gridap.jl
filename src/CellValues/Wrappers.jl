
struct CellValueFromArray{T,N,V<:AbstractArray{T,N}} <: IndexCellValue{T,N}
  v::V
end

@propagate_inbounds function getindex(self::CellValueFromArray,cell::Int)
  @inbounds self.v[cell]
end

@propagate_inbounds function getindex(self::CellValueFromArray{T,N},I::Vararg{Int,N}) where {T,N}
  @inbounds self.v[I...]
end

size(self::CellValueFromArray) = size(self.v)

IndexStyle(::Type{CellValueFromArray{T,N,V}}) where {T,N,V} = IndexStyle(V)

struct CellArrayFromArrayOfArrays{T,N,D,A,C<:AbstractArray{A,D}} <: IndexCellArray{T,N,A,D}
  c::C
end

function CellArrayFromArrayOfArrays(c::AbstractArray{A,D}) where {A<:AbstractArray,D}
  T = eltype(A)
  N = ndims(A)
  C = typeof(c)
  CellArrayFromArrayOfArrays{T,N,D,A,C}(c)
end

@propagate_inbounds function getindex(self::CellArrayFromArrayOfArrays,cell::Int)
  @inbounds self.c[cell]
end

@propagate_inbounds function getindex(self::CellArrayFromArrayOfArrays{T,N,D},I::Vararg{Int,D}) where {T,N,D}
  @inbounds self.c[I...]
end

size(self::CellArrayFromArrayOfArrays) = size(self.c)

IndexStyle(::Type{CellArrayFromArrayOfArrays{T,N,D,A,C}}) where {T,N,D,A,C} = IndexStyle(C)
