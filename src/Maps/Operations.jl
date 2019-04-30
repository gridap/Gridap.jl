# Unary operations

for op in (:+, :-)
  @eval begin
    function ($op)(a::Map)
      MapFromBroadcastUnaryOp($op,a)
    end
  end
end

"""
Abstract type representing a Map resulting from a generic unary operation
This includes, e.g., broadcasted operators, reductions, expansions, etc.
"""
abstract type MapFromUnaryOp{S,M,T,N} <: Map{S,M,T,N} end

inputmap(::MapFromUnaryOp) = @abstractmethod

computevals!(::MapFromUnaryOp, a, v) = @abstractmethod

computesize(::MapFromUnaryOp, sa) = @abstractmethod

cachedarray(::MapFromUnaryOp) = @abstractmethod

function return_size(
  self::MapFromUnaryOp{S,M}, s::NTuple{M,Int}) where {S,M}
  ma = inputmap(self)
  sa = return_size(ma,s)
  computesize(self,sa)
end

function evaluate!(
  self::MapFromUnaryOp{S,M,T,N},
  points::AbstractArray{S,M}, v::AbstractArray{T,N}) where {S,M,T,N}
  cache = cachedarray(self)
  ma = inputmap(self)
  sa = return_size(ma,size(points))
  setsize!(cache,sa)
  evaluate!(ma,points,cache)
  computevals!(self,cache,v)
end

gradient(::MapFromUnaryOp) = @notimplemented

"""
Map generated from a `Map` and a unary operator in broadcasted form
"""
struct MapFromBroadcastUnaryOp{
  O<:Function,S,M,T,N,U,A<:Map{S,M,U,N}} <: MapFromUnaryOp{S,M,T,N}
  op::O
  ma::A
  cache::CachedArray{U,N,Array{U,N}}
end

function MapFromBroadcastUnaryOp(
  op::Function, ma::Map{S,M,U,N}) where {S,M,U,N}
  cache = CachedArray(U,N)
  O = typeof(op)
  A = typeof(ma)
  T = Base._return_type(op,Tuple{U})
  MapFromBroadcastUnaryOp{O,S,M,T,N,U,A}(op,ma,cache)
end

inputmap(self::MapFromBroadcastUnaryOp) = self.ma

function computevals!(self::MapFromBroadcastUnaryOp, a, v)
  v .= self.op.(a)
end

computesize(self::MapFromBroadcastUnaryOp, sa) = sa

cachedarray(self::MapFromBroadcastUnaryOp) = self.cache

#"""
#Map generated from a `Map` and a unary operator
#"""
#struct MapFromUnaryOp{O,A,S,M,T,N} <: Map{S,M,T,N}
#  op::O
#  a::A
#end
#
#for op in (:+, :-)
#  @eval begin
#    function ($op)(a::Map)
#      MapFromUnaryOp($op,a)
#    end
#  end
#end
#
#function MapFromUnaryOp(op::Function,a::Map{S,M,T,N}) where {S,M,T,N}
#  O = typeof(op)
#  A = typeof(a)
#  R = Base._return_type(op,Tuple{T})
#  MapFromUnaryOp{O,A,S,M,R,N}(op,a)
#end
#
#function evaluate(self::MapFromUnaryOp{O,A,S,M,T,N},input::AbstractArray{S,M}) where {O,A,S,M,T,N}
#  avals = evaluate(self.a,input)
#  self.op(avals)
#end

# Binary operations

"""
Map generated from a binary operator over two `Map`
"""
struct MapFromBinaryOp{O,A,B,S,M,T,N} <: Map{S,M,T,N}
  op::O
  a::A
  b::B
end

for op in (:+, :-, :*, :/)
  @eval begin
    function ($op)(a::Map,b::Map)
      MapFromBinaryOp($op,a,b)
    end
  end
end

function inner(a::Field{D,T},b::Field{D,T}) where {D,T}
  MapFromBinaryOp(varinner,a,b)
end

function inner(a::Basis{D,T},b::Field{D,T}) where {D,T}
  MapFromBinaryOp(varinner,a,b)
end

function inner(a::Basis{D,T},b::Basis{D,T}) where {D,T}
  MapFromBinaryOp(varinner,a,b)
end

function expand(a::Basis,b::AbstractVector)
  println("No here")
  FieldFromExpand(a,b)
end


function MapFromBinaryOp(op::Function,a::Map{S,M,TA,NA},b::Map{S,M,TB,NB}) where {S,M,TA,NA,TB,NB}
  O = typeof(op)
  A = typeof(a)
  B = typeof(b)
  R = Base._return_type(op,Tuple{TA,TB})
  N = max(NA,NB)
  MapFromBinaryOp{O,A,B,S,M,R,N}(op,a,b)
end

function evaluate(self::MapFromBinaryOp{O,A,B,S,M,T,N},input::AbstractArray{S,M}) where {O,A,B,S,M,T,N}
  avals = evaluate(self.a,input)
  bvals = evaluate(self.b,input)
  self.op.(avals,bvals)
end

struct FieldFromExpand{D,S,R,T<:FieldValue} <: Field{D,T}
  basis::Basis{D,S}
  coeffs::AbstractVector{R}
end

function FieldFromExpand(basis::Basis{D,S},coeffs::AbstractVector{R}) where {D,S,R}
  T = Base._return_type(outer,Tuple{S,R})
  FieldFromExpand{D,S,R,T}(basis,coeffs)
end

function evaluate(self::FieldFromExpand{D},points::AbstractVector{Point{D}}) where D
  basisvals = evaluate(self.basis,points)
  expand(basisvals,self.coeffs)
end

function gradient(self::FieldFromExpand)
  gradbasis = gradient(self.basis)
  FieldFromExpand(gradbasis,self.coeffs)
end

varinner(a::FieldValue,b::FieldValue) = inner(a,b)

const BasisValues{T} = AbstractArray{T,2} where T <: FieldValue

varinner(a::BasisValues{T},b::T) where T<:FieldValue = inner(a,newaxis(b,dim=1))

varinner(a::BasisValues{T},b::BasisValues{T}) where T<:FieldValue = inner(newaxis(a,dim=2),newaxis(b,dim=1))

# expand(a::BasisValues,b) = cellsum(outer.(a,newaxis(b,dim=2)),dim=1)
function expand(a::AbstractArray{T,2}, b::AbstractVector{S}) where {T,S}
  num_dofs = length(b)
  @assert size(a)[1] == length(b)
  num_points = size(a)[end]
  R = Base._return_type(outer,Tuple{T,S})
  c = zeros(R,num_points)
  for i in 1:num_dofs
    for j in 1:num_points
      c[j] += a[i,j]*b[i]
    end
  end
  return c
end
# @santiagobadia : Think about efficiency improvements or other methods to
# implement this function

function newaxis(self::AbstractArray;dim::Int)
  s = [ v for v in size(self)]
  insert!(s,dim,1)
  shape = tuple(s...)
  return copy(reshape(self,shape))
end
