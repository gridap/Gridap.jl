
struct ConstantCellValue{T} <: IndexCellValue{T}
  value::T
  length::Int
end

size(self::ConstantCellValue) = (self.length,)

getindex(self::ConstantCellValue,cell::Int) = self.value

IndexStyle(::Type{ConstantCellValue{T}} where T) = IndexLinear()

const ConstantCellArray{T,N} = ConstantCellValue{<:AbstractArray{T,N}} where {T,N}

cellsize(self::ConstantCellArray) = size(self.value)
