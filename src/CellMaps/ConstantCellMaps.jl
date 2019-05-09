module ConstantCellMaps

using Numa.Helpers
using Numa.Maps
using Numa.CellValues
using Numa.CellMaps
using Numa.CellValues.ConstantCellValues

import Numa: evaluate
import Numa: gradient
export ConstantCellMap

const ConstantCellMap{S,M,T,N,R<:Map{S,M,T,N}} = ConstantCellValue{R}

function ConstantCellMap(a::Map{S,M,T,N},l::Int) where {S,M,T,N}
  R = typeof(a)
  ConstantCellMap{S,M,T,N,R}(a,l)
end

"""
Evaluate a `ConstantCellMap` on a set of points represented with a
`CellArray{S,M}`
"""
function evaluate(self::ConstantCellMap{S,M}, points::CellArray{S,M}) where {S,M}
  @notimplemented
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
