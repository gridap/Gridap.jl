
using Gridap.CellValues: IndexCellValue, IndexCellArray
import Base: size, getindex, IndexStyle, length
import Gridap.CellValues: cellsize
using StaticArrays

struct TestCellValue{T} <: IndexCellValue{T,1}
  a::T
  l::Int
end

size(self::TestCellValue) = (self.l,)

getindex(self::TestCellValue,cell::Int) = self.a

IndexStyle(::Type{TestCellValue{T}} where T) = IndexLinear()

struct TestCellArray{T,N} <: IndexCellArray{T,N,Array{T,N},1}
  a::Array{T,N}
  l::Int
end

size(self::TestCellArray) = (self.l,)

IndexStyle(::Type{TestCellArray{T,N}} where {T,N}) = IndexLinear()

getindex(self::TestCellArray,cell::Int) = self.a

cellsize(self::TestCellArray) = size(self.a)
