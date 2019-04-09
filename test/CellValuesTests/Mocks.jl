
using Numa.CellValues: IndexCellValue, IndexCellArray
import Base: size, getindex, IndexStyle, length
import Numa.CellValues: cellsize

struct TestCellValue{T} <: IndexCellValue{T,1}
  a::T
  l::Int
end

size(self::TestCellValue) = (self.l,)

getindex(self::TestCellValue,cell::Int) = self.a

IndexStyle(::Type{TestCellValue{T}} where T) = IndexLinear()


struct TestCellArray{T,N} <: IndexCellArray{T,N}
  a::Array{T,N}
  l::Int
end

length(self::TestCellArray) = self.l

getindex(self::TestCellArray,cell::Int) = self.a

cellsize(self::TestCellArray) = size(self.a)
