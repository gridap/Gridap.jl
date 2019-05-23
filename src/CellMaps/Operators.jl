module Operations

using Gridap
using Gridap.Helpers
using Gridap.Maps
using Gridap.CellValues
using Gridap.FieldValues
using Gridap.CellValues.Operations: CellValueFromUnaryOp
using Gridap.CellValues.Operations: CellValueFromBinaryOp
using Gridap.Maps: FieldFromExpand
using Gridap.CellMaps

import Gridap: evaluate, gradient
import Gridap: return_size
import Gridap.Maps: compose
import Gridap.Maps: lincomb
import Gridap.Maps: varinner
import Gridap.Maps: attachgeomap
import Base: iterate, length
import Gridap: apply

# Unary operations

# (TODO Remark) This implementation reuses the code from CellValueFromUnaryOp
# which is a nice thing, but, when iterating over CellMapFromUnaryOp, new 
# objects of type MapFromUnaryOp will create at each iteration ...
# The workaround would be to make MapFromUnaryOp mutable, and reimplement
# CellMapFromUnaryOp so that only a MapFromUnaryOp is created at the beginning
# and then at each iteration we reset the internal Map in MapFromUnaryOp

const CellMapFromUnaryOp{
  S,M,T,N,O,R,C} = CellValueFromUnaryOp{R,O,C} where {R<:Map{S,M,T,N},O<:Function}

function CellMapFromUnaryOp(op::Function,cellmap::CellMap)
  CellValueFromUnaryOp(op,cellmap)
end

function evaluate(self::CellMapFromUnaryOp{S,M,T,N},a::CellArray{S,M}) where {S,M,T,N}
  v = self.values
  x = evaluate(v,a)
  self.op(x)
end

function gradient(::CellMapFromUnaryOp{S,M,T,N}) where {S,M,T<:FieldValue,N}
  @notimplemented
end

function return_size(
  self::CellMapFromUnaryOp{S,M,T,N},s::NTuple{M,Int}) where {S,M,T,N}
  v, vnext = iterate(self)
  return_size(v,s)
end

# Binary operations

# (TODO remark) This implementation has the same pros and cons as CellMapFromUnaryOp

const CellMapFromBinaryOp{
  S,M,T,N,O,R,A,B} = CellValueFromBinaryOp{R,O,A,B} where {O<:Function,R<:Map{S,M,T,N}}

function CellMapFromBinaryOp(op::Function,a::CellMap{S,M},b::CellMap{S,M}) where {S,M}
  CellValueFromBinaryOp(op,a,b)
end

function evaluate(self::CellMapFromBinaryOp{S,M,T,N},q::CellArray{S,M}) where {S,M,T,N}
  a = self.a
  b = self.b
  xa = evaluate(a,q)
  xb = evaluate(b,q)
  self.op(xa,xb)
end

function gradient(::CellMapFromBinaryOp{S,M,T,N}) where {S,M,T<:FieldValue,N}
  @notimplemented
end

function gradient(v::CellValueFromBinaryOp{<:Map,typeof(+)})
  ga = gradient(v.a)
  gb = gradient(v.b)
  ga + gb
end

function gradient(v::CellValueFromBinaryOp{<:Map,typeof(-)})
  ga = gradient(v.a)
  gb = gradient(v.b)
  ga - gb
end

function return_size(
  self::CellMapFromBinaryOp{S,M,T,N},s::NTuple{M,Int}) where {S,M,T,N}
  v, vnext = iterate(self)
  return_size(v,s)
end

function varinner( a::CellField{D,T}, b::CellField{D,T}) where {D,T}
  _varinner(a,b)
end

function varinner( a::CellBasis{D,T}, b::CellField{D,T}) where {D,T}
  _varinner(a,b)
end

function varinner( a::CellBasis{D,T}, b::CellBasis{D,T}) where {D,T}
  _varinner(a,b)
end

function varinner(a::CellVector{T},b::CellVector{T}) where T
  inner(a,b)
end

function varinner(a::CellMatrix{T},b::CellVector{T}) where T
  inner(a,cellnewaxis(b,dim=1))
end

function varinner(a::CellMatrix{T},b::CellMatrix{T}) where T
  inner(cellnewaxis(a,dim=2),cellnewaxis(b,dim=1))
end

function _varinner(a,b)
  @assert length(a) == length(b)
  CellMapFromBinaryOp(varinner,a,b)
end

# lincomb

function lincomb(basis::CellBasis,coeffs::CellVector)
  CellFieldFromExpand(basis,coeffs)
end

struct CellFieldFromExpand{
  D,T,A<:CellBasis{D},B<:CellVector,R<:FieldFromExpand} <: IterCellField{D,T,R}
  basis::A
  coeffs::B
end

function CellFieldFromExpand(
  basis::CellBasis{D,S},coeffs::CellVector{C}) where {D,S,C}
  @assert length(basis) == length(coeffs)
  T = Base._return_type(outer,Tuple{S,C})
  A = typeof(basis)
  B = typeof(coeffs)
  R = FieldFromExpand{D,S,C,T,eltype(basis),eltype(coeffs)}
  CellFieldFromExpand{D,T,A,B,R}(basis,coeffs)
end

function length(this::CellFieldFromExpand)
  @assert length(this.basis) == length(this.coeffs)
  length(this.basis)
end

@inline function iterate(this::CellFieldFromExpand)
  bnext = iterate(this.basis)
  cnext = iterate(this.coeffs)
  iteratekernel(this,bnext,cnext)
end

@inline function iterate(this::CellFieldFromExpand,state)
  v, bstate, cstate = state
  bnext = iterate(this.basis,bstate)
  cnext = iterate(this.coeffs,cstate)
  iteratekernel(this,bnext,cnext)
end

function iteratekernel(this::CellFieldFromExpand,bnext,cnext)
  if bnext === nothing; return nothing end
  if cnext === nothing; return nothing end
  b, bstate = bnext
  c, cstate = cnext
  v = FieldFromExpand(b,c) #TODO This allocates a new object
  state = (v, bstate, cstate)
  (v, state)
end

function evaluate(self::CellFieldFromExpand{D},points::CellPoints{D}) where D
  basisvals = evaluate(self.basis,points)
  lincomb(basisvals,self.coeffs)
end

function lincomb(a::CellBasisValues,b::CellFieldValues)
  cellsum(outer(a,cellnewaxis(b,dim=2)),dim=1)
end

function gradient(self::CellFieldFromExpand)
  gradbasis = gradient(self.basis)
  CellFieldFromExpand(gradbasis,self.coeffs)
end

return_size(self::CellFieldFromExpand,s::Tuple{Int}) = s

end # module Operations


## Unary operations
#
#for op in (:+, :-)
#  @eval begin
#    function ($op)(a::CellMap)
#      CellMapFromUnaryOp($op,a)
#    end
#  end
#end
#
#struct CellMapFromUnaryOp{O,A,S,M,T,N} <: IterCellMap{S,M,T,N}
#  op::O
#  a::A
#end
#
#@inline function Base.iterate(this::CellMapFromUnaryOp{O,A,S,M,T,N}) where {O,A,S,M,T,N}
#  anext = iterate(this.a)
#  iteratekernel(this,anext)
#end
#
#@inline function Base.iterate(this::CellMapFromUnaryOp{O,A,S,M,T,N},state) where {O,A,S,M,T,N}
#  v, astate = state
#  anext = iterate(this.a,astate)
#  iteratekernel(this,anext)
#end
#
#function iteratekernel(this::CellMapFromUnaryOp,anext)
#  if anext === nothing; return nothing end
#  a, astate = anext
#  v = MapFromUnaryOp(this.op,a)
#  state = (v, astate)
#  (v, state)
#end
#
#function CellMapFromUnaryOp(op::Function,a::CellMap{S,M,T,N}) where {S,M,T,N}
#  O = typeof(op)
#  A = typeof(a)
#  R = Base._return_type(op,Tuple{T})
#  CellMapFromUnaryOp{O,A,S,M,R,N}(op,a)
#end
#
#function evaluate(self::CellMapFromUnaryOp{O,A,S,M,T,N},input::CellArray{S,M}) where {O,A,S,M,T,N}
#  avals = evaluate(self.a,input)
#  self.op(avals)
#end
#
## Binary operations
#
#for op in (:+, :-, :*, :/)
#  @eval begin
#    function ($op)(a::CellMap,b::CellMap)
#      CellMapFromBinaryOp($op,a,b)
#    end
#  end
#end
#
#function inner(a::CellField{D,T},b::CellField{D,T}) where {D,T}
#  CellMapFromBinaryOp(varinner,a,b)
#end
#
#function inner(a::CellBasis{D,T},b::CellField{D,T}) where {D,T}
#  CellMapFromBinaryOp(varinner,a,b)
#end
#
#function inner(a::CellBasis{D,T},b::CellBasis{D,T}) where {D,T}
#  CellMapFromBinaryOp(varinner,a,b)
#end
#
#function expand(a::CellBasis,b::CellVector)
#  CellFieldFromExpand(a,b)
#end
#
#struct CellMapFromBinaryOp{O,A,B,S,M,T,N} <: IterCellMap{S,M,T,N}
#  op::O
#  a::A
#  b::B
#end
#
#@inline function Base.iterate(this::CellMapFromBinaryOp{O,A,B,S,M,T,N}) where {O,A,B,S,M,T,N}
#  anext = iterate(this.a)
#  bnext = iterate(this.b)
#  iteratekernel(this,anext,bnext)
#end
#
#@inline function Base.iterate(this::CellMapFromBinaryOp{O,A,B,S,M,T,N},state) where {O,A,B,S,M,T,N}
#  v, astate, bstate = state
#  anext = iterate(this.a,astate)
#  bnext = iterate(this.b,bstate)
#  iteratekernel(this,anext,bnext)
#end
#
#function iteratekernel(this::CellMapFromBinaryOp,anext,bnext)
#  if anext === nothing; return nothing end
#  if bnext === nothing; return nothing end
#  a, astate = anext
#  b, bstate = bnext
#  v = MapFromBinaryOp(this.op,a,b)
#  state = (v, astate, bstate)
#  (v, state)
#end
#
#function CellMapFromBinaryOp(op::Function,a::CellMap{S,M,TA,NA},b::CellMap{S,M,TB,NB}) where {S,M,TA,NA,TB,NB}
#  O = typeof(op)
#  A = typeof(a)
#  B = typeof(b)
#  R = Base._return_type(op,Tuple{TA,TB})
#  N = max(NA,NB)
#  CellMapFromBinaryOp{O,A,B,S,M,R,N}(op,a,b)
#end
#
#function evaluate(self::CellMapFromBinaryOp{O,A,B,S,M,T,N},input::CellArray{S,M}) where {O,A,B,S,M,T,N}
#  avals = evaluate(self.a,input)
#  bvals = evaluate(self.b,input)
#  self.op(avals,bvals)
#end
#
#struct CellFieldFromExpand{D,S,R,T<:FieldValue} <: IterCellField{D,T}
#  basis::CellBasis{D,S}
#  coeffs::CellVector{R}
#end
#
#length(this::CellFieldFromExpand) = length(this.basis)
#
#@inline function Base.iterate(this::CellFieldFromExpand{D,S,R,T}) where {D,S,R,T}
#  bnext = iterate(this.basis)
#  cnext = iterate(this.coeffs)
#  iteratekernel(this,bnext,cnext)
#end
#
#@inline function Base.iterate(this::CellFieldFromExpand{D,S,R,T},state) where {D,S,R,T}
#  v, bstate, cstate = state
#  bnext = iterate(this.basis,bstate)
#  cnext = iterate(this.coeffs,cstate)
#  iteratekernel(this,bnext,cnext)
#end
#
#function iteratekernel(this::CellFieldFromExpand,bnext,cnext)
#  if bnext === nothing; return nothing end
#  if cnext === nothing; return nothing end
#  b, bstate = bnext
#  c, cstate = cnext
#  v = FieldFromExpand(b,c)
#  state = (v, bstate, cstate)
#  (v, state)
#end
#
#function CellFieldFromExpand(basis::CellBasis{D,S},coeffs::CellVector{R}) where {D,S,R}
#  T = Base._return_type(outer,Tuple{S,R})
#  CellFieldFromExpand{D,S,R,T}(basis,coeffs)
#end
#
#function evaluate(self::CellFieldFromExpand{D},points::CellPoints{D}) where D
#  basisvals = evaluate(self.basis,points)
#  expand(basisvals,self.coeffs)
#end
#
#function gradient(self::CellFieldFromExpand)
#  gradbasis = gradient(self.basis)
#  CellFieldFromExpand(gradbasis,self.coeffs)
#end
#
#varinner(a::FieldValue,b::FieldValue) = inner(a,b)
#
#varinner(a::CellFieldValues{T},b::CellFieldValues{T}) where T = inner(a,b)
#
#varinner(a::CellBasisValues{T},b::CellFieldValues{T}) where T = inner(a,cellnewaxis(b,dim=1))
#
#varinner(a::CellBasisValues{T},b::CellBasisValues{T}) where T = inner(cellnewaxis(a,dim=2),cellnewaxis(b,dim=1))
#
#expand(a::CellBasisValues,b::CellFieldValues) = cellsum(outer(a,cellnewaxis(b,dim=2)),dim=1)
#
## Composition
#
#function compose(f::Function,g::CellField{D,S}) where {D,S}
#  CellFieldFromCompose(f,g)
#end
#
#function compose(f::Function,g::CellGeomap{D,Z},u::CellField{Z,S}) where {D,Z,S}
#  CellFieldFromComposeExtended(f,g,u)
#end
#
#(∘)(f::Function,g::CellField) = compose(f,g)
#
#struct CellFieldFromCompose{D,O,C<:CellField{D},T} <: IterCellField{D,T}
#  a::C
#  op::O
#end
#
#length(this::CellFieldFromCompose) = length(this.a)
#
#@inline function Base.iterate(this::CellFieldFromCompose{D,O,C,T}) where {D,O,C,T}
#  anext = iterate(this.a)
#  iteratekernel(this,anext)
#end
#
#@inline function Base.iterate(this::CellFieldFromCompose{D,O,C,T},state) where {D,O,C,T}
#  v, astate = state
#  anext = iterate(this.a,astate)
#  iteratekernel(this,anext)
#end
#
#function iteratekernel(this::CellFieldFromCompose,anext)
#  if anext === nothing; return nothing end
#  a, astate = anext
#  v = FieldFromCompose(this.op,a)
#  state = (v, astate)
#  (v, state)
#end
#
#function CellFieldFromCompose(f::Function,g::CellField{D,S}) where {D,S}
#  O = typeof(f)
#  C = typeof(g)
#  T = Base._return_type(f,Tuple{S})
#  CellFieldFromCompose{D,O,C,T}(g,f)
#end
#
#function evaluate(self::CellFieldFromCompose{D},points::CellPoints{D}) where D
#  avals = evaluate(self.a,points)
#  composekernel(self,avals)
#end
#
#function composekernel(self::CellFieldFromCompose,avals::CellArray)
#  CellArrayFromBroadcastUnaryOp(self.op,avals)
#end
#
#function composekernel(self::CellFieldFromCompose,avals::ConstantCellArray)
#  b = broadcast(self.op,celldata(avals))
#  ConstantCellValue(b,length(avals))
#end
#
#function gradient(self::CellFieldFromCompose)
#  gradop = gradient(self.op)
#  CellFieldFromCompose(gradop,self.a)
#end
#
#struct CellFieldFromComposeExtended{D,O,G<:CellGeomap{D},U<:CellField,T} <: IterCellField{D,T}
#  f::O
#  g::G
#  u::U
#end
#
#length(this::CellFieldFromComposeExtended) = length(this.u)
#
#@inline function Base.iterate(this::CellFieldFromComposeExtended{D,O,G,U,T}) where {D,O,G,U,T}
#  gnext = iterate(this.g)
#  unext = iterate(this.u)
#  iteratekernel(this,gnext,unext)
#end
#
#@inline function Base.iterate(this::CellFieldFromComposeExtended{D,O,G,U,T},state) where {D,O,G,U,T}
#  v, gstate, ustate = state
#  gnext = iterate(this.g,astate)
#  unext = iterate(this.u,bstate)
#  iteratekernel(this,gnext,unext)
#end
#
#function iteratekernel(this::CellFieldFromComposeExtended,gnext,unext)
#  if gnext === nothing; return nothing end
#  if unext === nothing; return nothing end
#  g, gstate = gnext
#  u, ustate = unext
#  v = FieldFromComposeExtended(this.op,g,u)
#  state = (v, gstate, ustate)
#  (v, state)
#end
#
#function CellFieldFromComposeExtended(f::Function,g::CellGeomap{D,Z},u::CellField{Z,S}) where {D,Z,S}
#  O = typeof(f)
#  G = typeof(g)
#  U = typeof(u)
#  T = Base._return_type(f,Tuple{Point{Z},S})
#  CellFieldFromComposeExtended{D,O,G,U,T}(f,g,u)
#end
#
#function evaluate(self::CellFieldFromComposeExtended{D},points::CellPoints{D}) where D
#  gvals = evaluate(self.g,points)
#  uvals = evaluate(self.u,gvals)
#  CellArrayFromBroadcastBinaryOp(self.f,gvals,uvals)
#end
#
#function gradient(self::CellFieldFromComposeExtended)
#  gradf = gradient(self.f)
#  CellFieldFromComposeExtended(gradf,self.g,self.u)
#end

# @santiagobadia: With the following structure + unary/binary operators we can
# express CellFieldFromComposeExtended (after CellMapFromCompose)
# struct CellMapFromComposition{Q,O,A<:CellMap{S,M,Q,O},B<:CellMap{Q,O,T,N}) <: IterCellMap{S,M,T,N}
#   a::A
#   b::B
# end
