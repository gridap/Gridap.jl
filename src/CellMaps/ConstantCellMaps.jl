module ConstantCellMaps

using Numa.Helpers
using Numa.Maps
using Numa.CellValues
using Numa.CellMaps
using Numa.CellValues.ConstantCellValues
using Numa.CellMaps.CellMapValues

import Numa: evaluate
import Numa: gradient
import Numa: return_size
export ConstantCellMap

const ConstantCellMap{S,M,T,N,R<:Map{S,M,T,N}} = ConstantCellValue{R}

function ConstantCellMap(a::Map{S,M,T,N},l::Int) where {S,M,T,N}
  R = typeof(a)
  ConstantCellMap{S,M,T,N,R}(a,l)
end

function return_size(m::ConstantCellMap{S,M},s::NTuple{M,Int}) where {S,M}
  return_size(m.value,s)
end

"""
Evaluate a `ConstantCellMap` on a set of points represented with a
`CellArray{S,M}`
"""
function evaluate(self::ConstantCellMap{S,M}, points::CellArray{S,M}) where {S,M}
  CellMapValue(self,points)
end

function evaluate(self::ConstantCellMap{S,M}, points::ConstantCellArray{S,M}) where {S,M}
  @assert length(self) == length(points)
  l = self.length
  x = points.value
  m = self.value
  y = evaluate(m,x)
  ConstantCellArray(y,l)
end

"""
Computes the gradient of a `ConstantCellMap`
"""
function gradient(self::ConstantCellMap)
  gradfield = gradient(self.value)
  ConstantCellMap(gradfield,self.length)
end


end # module ConstantCellMaps
