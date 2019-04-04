
struct ConstantCellValue{T} <: IndexCellValue{T}
  value::T
  length::Int
end

size(self::ConstantCellValue) = (self.length,)

getindex(self::ConstantCellValue,cell::Int) = self.value

IndexStyle(::Type{ConstantCellValue{T}} where T) = IndexLinear()

const ConstantCellArray{T,N} = ConstantCellValue{<:AbstractArray{T,N}} where {T,N}

const ConstantCellVector{T,N} = ConstantCellArray{T,1} where T

cellsize(self::ConstantCellArray) = size(self.value)

ConstantCellArray(a::AbstractArray,l) = ConstantCellValue(a,l)

ConstantCellVector(a::AbstractVector,l) = ConstantCellValue(a,l)

function (==)(a::ConstantCellValue,b::ConstantCellValue)
  a.value != b.value && return false
  a.length != b.length && return false
  return true
end

for op in (:+,:-,:*,:/,:(inner),:(outer))

  @eval begin
    function ($op)(a::ConstantCellValue,b::ConstantCellValue)
      @assert length(a) == length(b)
      c = broadcast($op,a.value,b.value)
      ConstantCellValue(c,a.length)
    end
  end

end

for op in (:+,:-,:(det),:(inv))

  @eval begin
    function ($op)(a::ConstantCellValue)
      c = broadcast($op,a.value)
      ConstantCellValue(c,a.length)
    end
  end

end

