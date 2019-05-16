
using Gridap.CellMaps.CellMapValues: CellMapValue

import Gridap: gradient, evaluate, return_size

include("../CellValuesTests/Mocks.jl")

const TestCellMap{S,M,T,N,R} = TestCellValue{R} where R<:Map{S,M,T,N}

function TestCellMap(m::Map,l::Int)
  TestCellValue(m,l)
end

function evaluate(self::TestCellMap{S,M,T,N},a::CellArray{S,M}) where {S,M,T,N}
  CellMapValue(self,a)
end

function gradient(self::TestCellMap{S,M,T,N}) where {S,M,T,N}
  gradfield = gradient(self.a)
  TestCellMap(gradfield,self.l)
end

function return_size(
  self::TestCellMap{S,M},s::NTuple{M,Int}) where {S,M}
  return_size(self.a,s)
end
