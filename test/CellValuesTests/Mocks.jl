
using Numa.CellValues: IndexCellValue
import Base: size, getindex, IndexStyle
import Numa.CellValues: cellsize

struct TestCellValue{T} <: IndexCellValue{T}
  a::T
  l::Int
end

size(self::TestCellValue) = (self.l,)

getindex(self::TestCellValue,cell::Int) = self.a

IndexStyle(::Type{TestCellValue{T}} where T) = IndexLinear()

cellsize(self::TestCellValue{<:AbstractArray}) = size(self.a)
