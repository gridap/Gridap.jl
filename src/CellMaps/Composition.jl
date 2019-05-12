module Composition

using Numa
using Numa.Helpers
using Numa.Maps
using Numa.CellValues
using Numa.CellValues.ConstantCellValues
using Numa.CellValues.Operations: CellArrayFromBroadcastUnaryOp
using Numa.FieldValues
using Numa.CellMaps
using Numa.Maps: FieldFromCompose

import Numa: evaluate, gradient
import Numa: return_size
import Numa.Maps: compose
import Numa.Maps: attachgeomap
import Base: iterate, length


function compose(f::Function,g::CellField)
  CellFieldFromCompose(f,g)
end

function compose(f::Function,g::CellGeomap{D,Z},u::CellField{Z}) where {D,Z}
  CellFieldFromComposeExtended(f,g,u)
end

(∘)(f::Function,g::CellField) = compose(f,g)

# TODO this struct can be merged with CellMapFromBinaryOp

struct CellFieldFromCompose{
  D,T,O,C<:CellField{D},R<:FieldFromCompose} <: IterCellField{D,T,R}
  a::C
  op::O
end

function CellFieldFromCompose(f::Function,g::CellField{D,S}) where {D,S}
  O = typeof(f)
  C = typeof(g)
  T = Base._return_type(f,Tuple{S})
  R = FieldFromCompose{D,O,T,S,eltype(g)}
  CellFieldFromCompose{D,T,O,C,R}(g,f)
end

length(this::CellFieldFromCompose) = length(this.a)

@inline function Base.iterate(this::CellFieldFromCompose)
  anext = iterate(this.a)
  iteratekernel(this,anext)
end

@inline function Base.iterate(this::CellFieldFromCompose,state)
  v, astate = state
  anext = iterate(this.a,astate)
  iteratekernel(this,anext)
end

function iteratekernel(this::CellFieldFromCompose,anext)
  if anext === nothing; return nothing end
  a, astate = anext
  v = FieldFromCompose(this.op,a) # TODO this allocates a new object
  state = (v, astate)
  (v, state)
end

function evaluate(self::CellFieldFromCompose{D},points::CellPoints{D}) where D
  avals = evaluate(self.a,points)
  composekernel(self,avals)
end

function composekernel(self::CellFieldFromCompose,avals::CellArray)
  CellArrayFromBroadcastUnaryOp(self.op,avals)
end

function composekernel(self::CellFieldFromCompose,avals::ConstantCellArray)
  b = broadcast(self.op,celldata(avals))
  ConstantCellValue(b,length(avals))
end

function gradient(self::CellFieldFromCompose)
  gradop = gradient(self.op)
  CellFieldFromCompose(gradop,self.a)
end

return_size(self::CellFieldFromCompose,s::Tuple{Int}) = s

end # module Composition
