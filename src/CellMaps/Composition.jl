module Composition

using Gridap
using Gridap.Helpers
using Gridap.Maps
using Gridap.CellValues
using Gridap.CellValues.ConstantCellValues
using Gridap.CellValues.Operations: CellArrayFromBroadcastUnaryOp
using Gridap.CellValues.Operations: CellArrayFromBroadcastBinaryOp
using Gridap.FieldValues
using Gridap.CellMaps
using Gridap.Maps: FieldFromCompose
using Gridap.Maps: FieldFromComposeExtended

import Gridap: evaluate, gradient
import Gridap: return_size
import Gridap: compose
import Base: iterate, length, ∘

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
  ConstantCellArray(b,length(avals))
end

function gradient(self::CellFieldFromCompose)
  gradop = gradient(self.op)
  CellFieldFromCompose(gradop,self.a)
end

return_size(self::CellFieldFromCompose,s::Tuple{Int}) = s

struct CellFieldFromComposeExtended{
  D,T,O,G<:CellGeomap{D},U<:CellField,R<:FieldFromComposeExtended} <: IterCellField{D,T,R}
  f::O
  g::G
  u::U
end

function CellFieldFromComposeExtended(
  f::Function,g::CellGeomap{D,Z},u::CellField{Z,S}) where {D,Z,S}
  @assert length(g) == length(u)
  O = typeof(f)
  G = typeof(g)
  U = typeof(u)
  T = Base._return_type(f,Tuple{Point{Z},S})
  R = FieldFromComposeExtended{D,T,O,Z,S,eltype(g),eltype(u)}
  CellFieldFromComposeExtended{D,T,O,G,U,R}(f,g,u)
end

length(this::CellFieldFromComposeExtended) = length(this.u)

@inline function Base.iterate(this::CellFieldFromComposeExtended)
  gnext = iterate(this.g)
  unext = iterate(this.u)
  iteratekernel(this,gnext,unext)
end

@inline function Base.iterate(this::CellFieldFromComposeExtended,state)
  v, gstate, ustate = state
  gnext = iterate(this.g,gstate)
  unext = iterate(this.u,ustate)
  iteratekernel(this,gnext,unext)
end

function iteratekernel(this::CellFieldFromComposeExtended,gnext,unext)
  if gnext === nothing; return nothing end
  if unext === nothing; return nothing end
  g, gstate = gnext
  u, ustate = unext
  v = FieldFromComposeExtended(this.f,g,u) # TODO this allocates a new object
  state = (v, gstate, ustate)
  (v, state)
end

function evaluate(self::CellFieldFromComposeExtended{D},points::CellPoints{D}) where D
  gvals = evaluate(self.g,points)
  uvals = evaluate(self.u,points)
  composekernel(self,gvals,uvals)
end

function composekernel(
  self::CellFieldFromComposeExtended,gvals::CellArray,uvals::CellArray)
  CellArrayFromBroadcastBinaryOp(self.f,gvals,uvals)
end

function composekernel(
  self::CellFieldFromComposeExtended,
  gvals::ConstantCellArray,
  uvals::ConstantCellArray)
  b = broadcast(self.f,celldata(gvals),celldata(uvals))
  ConstantCellArray(b,length(gvals))
end

function gradient(self::CellFieldFromComposeExtended)
  gradf = gradient(self.f)
  CellFieldFromComposeExtended(gradf,self.g,self.u)
end

return_size(self::CellFieldFromComposeExtended,s::Tuple{Int}) = s

# @santiagobadia: With the following structure + unary/binary operators we can
# express CellFieldFromComposeExtended (after CellMapFromCompose)
# struct CellMapFromComposition{Q,O,A<:CellMap{S,M,Q,O},B<:CellMap{Q,O,T,N}) <: IterCellMap{S,M,T,N}
#   a::A
#   b::B
# end

end # module Composition
