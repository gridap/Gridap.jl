module ConstantCellMaps

using Gridap.Helpers
using Gridap.FieldValues
using Gridap.Maps
using Gridap.CellValues
using Gridap.CellMaps
using Gridap.CellValues.ConstantCellValues
using Gridap.CellMaps.CellMapValues

import Gridap: evaluate
import Gridap: gradient
import Gridap: return_size
import Gridap.Maps: compose
import Gridap.Maps: lincomb
import Gridap.Maps: varinner
import Gridap.Maps: attachgeomap
export ConstantCellMap
export ConstantCellField
export ConstantCellBasis

const ConstantCellMap{S,M,T,N,R<:Map{S,M,T,N}} = ConstantCellValue{R}

function ConstantCellMap(a::Map{S,M,T,N},l::Int) where {S,M,T,N}
  R = typeof(a)
  ConstantCellMap{S,M,T,N,R}(a,l)
end

const ConstantCellField{D,T,R<:Field{D,T}} = ConstantCellMap{Point{D},1,T,1,R} where T<:FieldValue

const ConstantCellBasis{D,T,R<:Basis{D,T}} = ConstantCellMap{Point{D},1,T,2,R} where T<:FieldValue

function ConstantCellField(a::Field,l::Int)
  ConstantCellMap(a,l)
end

function ConstantCellBasis(a::Basis,l::Int)
  ConstantCellMap(a,l)
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

# Optimized operations

function compose(f::Function,g::ConstantCellField)
  m = compose(f,g.value)
  ConstantCellMap(m,g.length)
end

function compose(
  f::Function,
  g::ConstantCellField{D,Point{Z}},
  u::ConstantCellField{Z,T}) where {D,Z,T}

  @assert g.length == u.length
  m = compose(f,g.value,u.value)
  ConstantCellMap(m,g.length)
end

function varinner(
  a::ConstantCellField{D,T},
  b::ConstantCellField{D,T}) where {D,T}
  _varinner(a,b)
end

function varinner(
  a::ConstantCellBasis{D,T},
  b::ConstantCellField{D,T}) where {D,T}
  _varinner(a,b)
end

function varinner(
  a::ConstantCellBasis{D,T},
  b::ConstantCellBasis{D,T}) where {D,T}
  _varinner(a,b)
end

function _varinner(a,b)
  @assert a.length == b.length
  m = varinner(a.value,b.value)
  ConstantCellMap(m,a.length)
end

function lincomb(a::ConstantCellBasis,b::ConstantCellVector)
  @assert a.length == b.length
  m = lincomb(a.value,b.value)
  ConstantCellMap(m,a.length)
end

function attachgeomap(
  a::ConstantCellBasis{D}, b::ConstantCellField{D,Point{D}}) where D
  @assert a.length == b.length
  m = attachgeomap(a.value,b.value)
  ConstantCellMap(m,a.length)
end

end # module ConstantCellMaps
