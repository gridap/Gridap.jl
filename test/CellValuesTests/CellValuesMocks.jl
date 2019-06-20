module CellValuesMocks

using Gridap

export TestIterCellValue
export TestIndexCellValue
import Base: iterate
import Base: length
import Base: size
import Base: getindex
import Gridap: evaluate

struct TestIndexCellValue{T} <: IndexCellValue{T,1}
  a::T
  l::Int
end

size(self::TestIndexCellValue) = (self.l,)

getindex(self::TestIndexCellValue,cell::Int) = self.a

struct TestIterCellValue{T} <: IterCellValue{T}
  a::T
  l::Int
end

length(v::TestIterCellValue) = v.l

@inline function iterate(v::TestIterCellValue)
  i = 0
  if v.l <= i; return nothing; end
  (v.a,i+1)
end

@inline function iterate(v::TestIterCellValue,state)
  i = state
  if v.l <= i; return nothing; end
  (v.a,i+1)
end

end # module
