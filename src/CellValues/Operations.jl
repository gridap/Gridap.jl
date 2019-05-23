module Operations

using Gridap.Helpers
using Gridap.CellValues

using StaticArrays
using Base: @propagate_inbounds
using Base.Cartesian: @nloops, @nexprs, @nref
using Gridap.CachedArrays
using Gridap.Maps: newaxis_kernel!, newaxis_size

export cellsum
export cellnewaxis
export cellmean
import Base: +, -, *, /
import Base: ==
import LinearAlgebra: inv, det
import Gridap.FieldValues: inner, outer, meas

import Base: iterate
import Base: length
import Base: eltype
import Base: size
import Base: getindex
import Base: IndexStyle
import Gridap.CellValues: cellsize
import Gridap: apply

# Unary operations on CellValue

for op in (:+, :-, :(inv), :(det), :(meas))
  @eval begin
    function ($op)(a::CellValue)
      CellValueFromUnaryOp($op,a)
    end
    function ($op)(a::CellValue{<:SArray})
      CellValueFromUnaryOp($op,a)
    end
  end
end

function apply(op::Function,a::CellValue)
  CellValueFromUnaryOp(op,a)
end

function apply(op::Function,a::CellValue{<:SArray})
  CellValueFromUnaryOp(op,a)
end

struct CellValueFromUnaryOp{T,O<:Function,C<:CellValue} <: IterCellValue{T}
  op::O
  values::C
end

function CellValueFromUnaryOp(op::Function,values::CellValue{T}) where T
  O = typeof(op)
  C = typeof(values)
  S = Base._return_type(op,Tuple{T})
  CellValueFromUnaryOp{S,O,C}(op,values)
end

length(self::CellValueFromUnaryOp) = length(self.values)

@inline function iterate(self::CellValueFromUnaryOp)
  next = iterate(self.values)
  if next === nothing; return nothing end
  a, state = next
  (self.op(a), state)
end

@inline function iterate(self::CellValueFromUnaryOp,state)
  next = iterate(self.values,state)
  if next === nothing; return nothing end
  a, state = next
  (self.op(a), state)
end

# Binary operations on CellValue

function (==)(a::CellValue{T},b::CellValue{T}) where T
  _cell_value_eq_kernel(a,b)
end

function (==)(a::CellValue{T},b::CellValue{T}) where T<:SArray
  _cell_value_eq_kernel(a,b)
end

function _cell_value_eq_kernel(a,b)
  length(a) != length(b) && return false
  for (ai,bi) in zip(a,b)
    ai != bi && return false
  end
  return true
end

for op in (:+, :-, :*, :/, :(outer), :(inner))
  @eval begin
    function ($op)(a::CellValue,b::CellValue)
      CellValueFromBinaryOp($op,a,b)
    end
    function ($op)(a::CellValue{<:SArray},b::CellValue{<:SArray})
      CellValueFromBinaryOp($op,a,b)
    end
  end
end

function apply(op::Function,a::CellValue,b::CellValue)
  CellValueFromBinaryOp(op,a,b)
end

# Ancillary types

struct CellValueFromBinaryOp{T,O<:Function,A<:CellValue,B<:CellValue} <: IterCellValue{T}
  op::O
  a::A
  b::B
end

function CellValueFromBinaryOp(op::Function,a::CellValue{T},b::CellValue{S}) where {T,S}
  @assert length(a) == length(b)
  O = typeof(op)
  A = typeof(a)
  B = typeof(b)
  R = Base._return_type(op,Tuple{T,S})
  CellValueFromBinaryOp{R,O,A,B}(op,a,b)
end

function length(self::CellValueFromBinaryOp)
  @assert length(self.a) == length(self.b)
  length(self.a)
end

@inline function iterate(self::CellValueFromBinaryOp)
  anext = iterate(self.a)
  bnext = iterate(self.b)
  if anext === nothing; return nothing end
  if bnext === nothing; return nothing end
  a, astate = anext
  b, bstate = bnext
  state = (astate,bstate)
  (self.op(a,b), state)
end

@inline function iterate(self::CellValueFromBinaryOp,state)
  astate, bstate = state
  anext = iterate(self.a,astate)
  bnext = iterate(self.b,bstate)
  if anext === nothing; return nothing end
  if bnext === nothing; return nothing end
  a, astate = anext
  b, bstate = bnext
  state = (astate,bstate)
  (self.op(a,b), state)
end

# Unary operations on CellArray

for op in (:+, :-, :(inv), :(det), :(meas))
  @eval begin
    function ($op)(a::CellArray)
      CellArrayFromBroadcastUnaryOp($op,a)
    end
  end
end

function apply(op::Function,a::CellArray)
  CellArrayFromBroadcastUnaryOp(op,a)
end

function cellsum(self::CellArray{T,N};dim::Int) where {T,N}
  CellArrayFromCellSum{dim,N-1,typeof(self),T}(self)
end

function cellsum(self::CellArray{T,1};dim::Int) where T
  CellValueFromCellArrayReduce(sum,self)
end

function cellnewaxis(self::CellArray{T,N};dim::Int) where {T,N}
  CellArrayFromCellNewAxis{dim,typeof(self),T,N+1}(self)
end

mean(a) = sum(a)/length(a)

function cellmean(self::CellArray)
  CellValueFromCellArrayReduce(mean,self)
end

# Ancillary types

abstract type CellArrayFromUnaryOp{
  C<:CellArray,T,N} <: IterCellArray{T,N,CachedArray{T,N,Array{T,N}}} end

function inputcellarray(::CellArrayFromUnaryOp{C,T,N})::C  where {C,T,N}
  @abstractmethod
end

function computesize(::CellArrayFromUnaryOp{C,T,N}, asize) where {C,T,N}
  @abstractmethod
end

computevals!(::CellArrayFromUnaryOp, a, v) = @abstractmethod

Base.length(self::CellArrayFromUnaryOp) = length(inputcellarray(self))

cellsize(self::CellArrayFromUnaryOp) = computesize(self,cellsize(inputcellarray(self)))

@inline function Base.iterate(self::CellArrayFromUnaryOp{C,T,N}) where {C,T,N}
  u = Array{T,N}(undef,cellsize(self))
  v = CachedArray(u)
  anext = iterate(inputcellarray(self))
  if anext === nothing; return nothing end
  iteratekernel(self,anext,v)
end

@inline function Base.iterate(self::CellArrayFromUnaryOp,state)
  v, astate = state
  anext = iterate(inputcellarray(self),astate)
  if anext === nothing; return nothing end
  iteratekernel(self,anext,v)
end

function iteratekernel(self::CellArrayFromUnaryOp,anext,v)
  a, astate = anext
  vsize = computesize(self,size(a))
  setsize!(v,vsize)
  computevals!(self,a,v)
  state = (v, astate)
  (v,state)
end

struct CellArrayFromBroadcastUnaryOp{
  O<:Function,C<:CellArray,T,N} <: CellArrayFromUnaryOp{C,T,N}
  op::O
  a::C
end

function CellArrayFromBroadcastUnaryOp(op::Function,a::CellArray{T,N}) where {T,N}
  O = typeof(op)
  C = typeof(a)
  S = Base._return_type(op,Tuple{T})
  CellArrayFromBroadcastUnaryOp{O,C,S,N}(op,a)
end

computesize(::CellArrayFromBroadcastUnaryOp, asize) = asize

inputcellarray(self::CellArrayFromBroadcastUnaryOp) = self.a

function computevals!(self::CellArrayFromBroadcastUnaryOp, a, v)
  broadcast!(self.op,v,a)
end

struct CellArrayFromCellSum{A,N,C,T} <: CellArrayFromUnaryOp{C,T,N}
  a::C
end

inputcellarray(self::CellArrayFromCellSum) = self.a

function computesize(self::CellArrayFromCellSum{A},asize) where A
  cellsumsize(asize,Val(A))
end

function computevals!(::CellArrayFromCellSum{A}, a, v) where A
  cellsumvals!(a,v,Val(A))
end

struct CellArrayFromCellNewAxis{A,C,T,N} <: CellArrayFromUnaryOp{C,T,N}
  a::C
end

inputcellarray(self::CellArrayFromCellNewAxis) = self.a

function computesize(self::CellArrayFromCellNewAxis{A},asize::NTuple{M,Int}) where {A,M}
  newaxis_size(Val(A),asize)
end

function computevals!(::CellArrayFromCellNewAxis{D}, A::AbstractArray{T,M}, B::AbstractArray{T,N}) where {D,T,M,N}
  newaxis_kernel!(Val(D),A,B)
end

struct CellValueFromCellArrayReduce{T,O<:Function,C<:CellArray} <: IterCellValue{T}
  op::O
  values::C
end

function CellValueFromCellArrayReduce(op::Function,values::CellArray{T,N}) where {T,N}
  O = typeof(op)
  C = typeof(values)
  S = Base._return_type(op,Tuple{Array{T,N}})
  CellValueFromCellArrayReduce{S,O,C}(op,values)
end

length(self::CellValueFromCellArrayReduce) = length(self.values)

@inline function iterate(self::CellValueFromCellArrayReduce)
  next = iterate(self.values)
  if next === nothing; return nothing end
  a, state = next
  (self.op(a), state)
end

@inline function iterate(self::CellValueFromCellArrayReduce,state)
  next = iterate(self.values,state)
  if next === nothing; return nothing end
  a, state = next
  (self.op(a), state)
end

# Binary operations on CellArray

function (==)(a::CellArray{T,N},b::CellArray{T,N}) where {T,N}
  length(a) != length(b) && return false
  cellsize(a) != cellsize(b) && return false
  for (ai,bi) in zip(a,b)
    ai != bi && return false
  end
  return true
end

for op in (:+, :-, :*, :/, :(inner), :(outer))
  @eval begin
    function ($op)(a::CellArray,b::CellArray)
      CellArrayFromBroadcastBinaryOp($op,a,b)
    end
    function ($op)(a::CellArray,b::CellValue)
      CellArrayFromBroadcastBinaryOp($op,a,b)
    end
    function ($op)(a::CellValue,b::CellArray)
      CellArrayFromBroadcastBinaryOp($op,a,b)
    end
  end
end

# Ancillary types

abstract type CellArrayFromBinaryOp{
  A<:CellValue,B<:CellValue,T,N} <: IterCellArray{T,N,CachedArray{T,N,Array{T,N}}} end

function leftcellarray(::CellArrayFromBinaryOp{A,B,T,N})::A where {A,B,T,N}
  @abstractmethod
end

function rightcellarray(::CellArrayFromBinaryOp{A,B,T,N})::B where {A,B,T,N}
  @abstractmethod
end

computesize(::CellArrayFromBinaryOp, asize, bsize) = @abstractmethod

computevals!(::CellArrayFromBinaryOp, a, b, v) = @abstractmethod

function Base.length(self::CellArrayFromBinaryOp)
  @assert length(rightcellarray(self)) == length(leftcellarray(self))
  length(rightcellarray(self))
end

cellsize(self::CellArrayFromBinaryOp) = computesize(self,cellsize(leftcellarray(self)),cellsize(rightcellarray(self)))

@inline function Base.iterate(self::CellArrayFromBinaryOp{A,B,T,N}) where {A,B,T,N}
  u = Array{T,N}(undef,cellsize(self))
  v = CachedArray(u)
  anext = iterate(leftcellarray(self))
  if anext === nothing; return nothing end
  bnext = iterate(rightcellarray(self))
  if bnext === nothing; return nothing end
  iteratekernel(self,anext,bnext,v)
end

@inline function Base.iterate(self::CellArrayFromBinaryOp,state)
  v, astate, bstate = state
  anext = iterate(leftcellarray(self),astate)
  if anext === nothing; return nothing end
  bnext = iterate(rightcellarray(self),bstate)
  if bnext === nothing; return nothing end
  iteratekernel(self,anext,bnext,v)
end

function iteratekernel(self::CellArrayFromBinaryOp,anext,bnext,v)
  a, astate = anext
  b, bstate = bnext
  vsize = computesize(self,_custom_size(a),_custom_size(b))
  setsize!(v,vsize)
  computevals!(self,a,b,v)
  state = (v, astate, bstate)
  (v,state)
end

_custom_size(a) = size(a)

_custom_size(a::SArray) = ()

struct CellArrayFromBroadcastBinaryOp{O<:Function,T,N,A,B} <: CellArrayFromBinaryOp{A,B,T,N}
  op::O
  a::A
  b::B
end

function CellArrayFromBroadcastBinaryOp(op::Function,a::CellArray{T,N},b::CellArray{S,M}) where {T,S,N,M}
  O = typeof(op)
  A = typeof(a)
  B = typeof(b)
  R = Base._return_type(op,Tuple{T,S})
  L = max(N,M)
  CellArrayFromBroadcastBinaryOp{O,R,L,A,B}(op,a,b)
end

function CellArrayFromBroadcastBinaryOp(op::Function,a::CellArray{T,N},b::CellValue{S}) where {T,S,N}
  O = typeof(op)
  A = typeof(a)
  B = typeof(b)
  R = Base._return_type(op,Tuple{T,S})
  CellArrayFromBroadcastBinaryOp{O,R,N,A,B}(op,a,b)
end

function CellArrayFromBroadcastBinaryOp(op::Function,a::CellValue{T},b::CellArray{S,N}) where {T,S,N}
  O = typeof(op)
  A = typeof(a)
  B = typeof(b)
  R = Base._return_type(op,Tuple{T,S})
  CellArrayFromBroadcastBinaryOp{O,R,N,A,B}(op,a,b)
end

function CellArrayFromBroadcastBinaryOp(op::Function,a::CellArray{T,N},b::CellValue{S}) where {T,S<:SArray,N}
  O = typeof(op)
  A = typeof(a)
  B = typeof(b)
  R = Base._return_type(op,Tuple{T,S})
  CellArrayFromBroadcastBinaryOp{O,R,N,A,B}(op,a,b)
end

function CellArrayFromBroadcastBinaryOp(op::Function,a::CellValue{T},b::CellArray{S,N}) where {T<:SArray,S,N}
  O = typeof(op)
  A = typeof(a)
  B = typeof(b)
  R = Base._return_type(op,Tuple{T,S})
  CellArrayFromBroadcastBinaryOp{O,R,N,A,B}(op,a,b)
end

leftcellarray(self::CellArrayFromBroadcastBinaryOp) = self.a

rightcellarray(self::CellArrayFromBroadcastBinaryOp) = self.b

function computesize(::CellArrayFromBroadcastBinaryOp, asize, bsize)
  Base.Broadcast.broadcast_shape(asize,bsize)
end

function computevals!(self::CellArrayFromBroadcastBinaryOp, a, b, v)
  _custom_broadcast!(self.op,v,a,b)
end

function _custom_broadcast!(op,v,a,b)
  broadcast!(op,v,a,b)
end

function _custom_broadcast!(op,v::AbstractArray,a::SArray,b::AbstractArray)
  @inbounds for i in eachindex(b)
    v[i] = op(a,b[i])
  end
end

function _custom_broadcast!(op,v::AbstractArray,a::AbstractArray,b::SArray)
  @inbounds for i in eachindex(a)
    v[i] = op(a[i],b)
  end
end

function _custom_broadcast(op,a,b)
  broadcast(op,a,b)
end

function _custom_broadcast(
  op, a::SArray{Size,T}, b::AbstractArray{S,M}) where {Size,T,S,M}
  R = Base._return_type(op,Tuple{typeof(a),S})
  v = Array{R,M}(undef,size(b))
  _custom_broadcast!(op,v,a,b)
  v
end

function _custom_broadcast(
  op, a::AbstractArray{S,M}, b::SArray{Size,T}) where {Size,T,S,M}
  R = Base._return_type(op,Tuple{S,typeof(b)})
  v = Array{R,M}(undef,size(a))
  _custom_broadcast!(op,v,a,b)
  v
end

@generated function cellsumsize(asize::NTuple{N,Int},::Val{D}) where {N,D}
  @assert N > 0
  @assert D <= N
  str = join([ "asize[$i]," for i in 1:N if i !=D ])
  Meta.parse("($str)")
end

@generated function cellsumvals!(A::AbstractArray{T,N},B,::Val{D}) where {T,N,D}
  @assert N > 0
  @assert D <= N
  quote
    @nloops $(N-1) b B begin
      (@nref $(N-1) B b) = zero(T)
    end
    @nloops $N a A begin
      @nexprs $(N-1) j->(b_j = a_{ j < $D ? j : j+1  } )
      (@nref $(N-1) B b) += @nref $N A a 
    end
  end    
end

end # module Operations
