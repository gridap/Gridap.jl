module Gallery

using Gridap
using Gridap.Helpers
using Base: @propagate_inbounds, @pure

export CellValueFromArray
import Base: iterate
import Base: length
import Base: eltype
import Base: size
import Base: getindex
import Base: IndexStyle

struct CellValueFromArray{T,N,V<:AbstractArray{T,N}} <: IndexCellValue{T,N}
  v::V
end

@propagate_inbounds function getindex(
  self::CellValueFromArray{T,N},I::Vararg{Int,N}) where {T,N}
  @inbounds self.v[I...]
end

size(self::CellValueFromArray) = size(self.v)

@pure IndexStyle(::Type{CellValueFromArray{T,N,V}}) where {T,N,V} = IndexStyle(V)

end # module Gallery
