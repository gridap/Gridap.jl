# Unary operations

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

for op in (:+, :-)
  @eval begin
    function ($op)(a::Map)
      MapFromBroadcastUnaryOp($op,a)
    end
  end
end

"""
Map holding the result of applying the newaxis operation to another Map
"""
struct MapFromNewAxis{D,S,M,T,N,K,A<:Map{S,M,T,K}} <: MapFromUnaryOp{S,M,T,N}
  ma::A
  cache::CachedArray{T,K,Array{T,K}}
end

function MapFromNewAxis(m::Map{S,M,T,K}, dim::Integer) where {S,M,T,K}
  D = dim
  N = K+1
  A = typeof(m)
  cache = CachedArray(T,K)
  MapFromNewAxis{D,S,M,T,N,K,A}(m,cache)
end

inputmap(m::MapFromNewAxis) = m.ma

function computevals!(::MapFromNewAxis{D}, a, v) where D
  newaxis_kernel!(Val(D),a,v)
end

computesize(::MapFromNewAxis{D}, sa) where {D} = newaxis_size(Val(D),sa)

cachedarray(m::MapFromNewAxis) = m.cache

function newaxis(m::Map;dim::Integer)
  MapFromNewAxis(m,dim)
end

@generated function newaxis_size(::Val{A},asize::NTuple{M,Int}) where {A,M}
  @assert A <= M+1
  str = ["asize[$i]," for i in 1:M]
  insert!(str,A,"1,")
  Meta.parse("($(join(str)))")
end

@generated function newaxis_kernel!(::Val{D}, A::AbstractArray{T,M}, B::AbstractArray{T,N}) where {D,T,M,N}
  @assert N == M + 1
  @assert D <= N
  quote
    @nloops $M a A begin
      @nexprs $N j->(b_j = j == $D ? 1 : a_{ j < $D ? j : j-1 } )
      (@nref $N B b) = @nref $M A a
    end
  end
end

function lincomb(basis::Basis{D,S},coeffs::AbstractVector{R}) where {D,S,R}
  FieldFromExpand(basis,coeffs)
end

"""
Field resulting from the linear combination of a Basis and a coefficient vector
"""
struct FieldFromExpand{
  D,S,R,T,B<:Basis{D,S},C<:AbstractVector{R}} <: MapFromUnaryOp{Point{D},1,T,1}
  basis::B
  coeffs::C
  cache::CachedArray{S,2,Array{S,2}}
end

function FieldFromExpand(basis::Basis{D,S},coeffs::AbstractVector{R}) where {D,S,R}
  T = Base._return_type(outer,Tuple{S,R})
  B = typeof(basis)
  C = typeof(coeffs)
  cache = CachedArray(S,2)
  FieldFromExpand{D,S,R,T,B,C}(basis,coeffs,cache)
end

inputmap(m::FieldFromExpand) = m.basis

function computevals!(m::FieldFromExpand, a, v)
  b = m.coeffs
  lin_comb_kernel!(a,b,v)
end

function lin_comb_kernel!(a,b,v::AbstractArray{T}) where T
  ndofs, npoints = size(a)
  @assert length(b) == ndofs
  @assert length(v) == npoints
  for i in eachindex(v)
    v[i] = zero(T)
  end
  for j in 1:npoints
    for i in 1:ndofs
      v[j] += outer(a[i,j],b[i])
    end
  end
end

computesize(::FieldFromExpand, sa) = (sa[2],)

cachedarray(m::FieldFromExpand) = m.cache

function gradient(self::FieldFromExpand)
  gradbasis = gradient(self.basis)
  FieldFromExpand(gradbasis,self.coeffs)
end

# Binary operations

"""
Abstract type representing a Map resulting from a generic binary operation
"""
abstract type MapFromBinaryOp{S,M,T,N} <: Map{S,M,T,N} end

leftmap(::MapFromBinaryOp) = @abstractmethod

rightmap(::MapFromBinaryOp) = @abstractmethod

computevals!(::MapFromBinaryOp, a, b, v) = @abstractmethod

computesize(::MapFromBinaryOp, sa, sb) = @abstractmethod

leftcachedarray(::MapFromBinaryOp) = @abstractmethod

rightcachedarray(::MapFromBinaryOp) = @abstractmethod

function return_size(
  self::MapFromBinaryOp{S,M}, s::NTuple{M,Int}) where {S,M}
  ma = leftmap(self)
  mb = rightmap(self)
  sa = return_size(ma,s)
  sb = return_size(mb,s)
  computesize(self,sa,sb)
end

function evaluate!(
  self::MapFromBinaryOp{S,M,T,N},
  points::AbstractArray{S,M}, v::AbstractArray{T,N}) where {S,M,T,N}
  acache = leftcachedarray(self)
  bcache = rightcachedarray(self)
  ma = leftmap(self)
  mb = rightmap(self)
  s = size(points)
  sa = return_size(ma,s)
  sb = return_size(mb,s)
  setsize!(acache,sa)
  setsize!(bcache,sb)
  evaluate!(ma,points,acache)
  evaluate!(mb,points,bcache)
  computevals!(self,acache,bcache,v)
end

gradient(::MapFromBinaryOp) = @notimplemented

"""
Map generated from two `Map` objects and a binary operator in broadcasted form
"""
struct MapFromBroadcastBinaryOp{
  O<:Function,S,M,T,N,U,I,V,J,
  A<:Map{S,M,U,I},B<:Map{S,M,V,J}} <: MapFromBinaryOp{S,M,T,N}
  op::O
  ma::A
  mb::B
  acache::CachedArray{U,I,Array{U,I}}
  bcache::CachedArray{V,J,Array{V,J}}
end

function MapFromBroadcastBinaryOp(
  op::Function, ma::Map{S,M,U,I}, mb::Map{S,M,V,J}) where {S,M,U,I,V,J}
  acache = CachedArray(U,I)
  bcache = CachedArray(V,J)
  O = typeof(op)
  A = typeof(ma)
  B = typeof(mb)
  T = Base._return_type(op,Tuple{U,V})
  N = max(I,J)
  MapFromBroadcastBinaryOp{O,S,M,T,N,U,I,V,J,A,B}(op,ma,mb,acache,bcache)
end

leftmap(self::MapFromBroadcastBinaryOp) = self.ma

rightmap(self::MapFromBroadcastBinaryOp) = self.mb

function computevals!(self::MapFromBroadcastBinaryOp, a, b, v)
  v .= self.op.(a,b)
end

function computesize(self::MapFromBroadcastBinaryOp, sa, sb)
  Base.Broadcast.broadcast_shape(sa,sb)
end

leftcachedarray(self::MapFromBroadcastBinaryOp) = self.acache

rightcachedarray(self::MapFromBroadcastBinaryOp) = self.bcache

for op in (:+, :-, :*, :/, :(inner), :(outer))
  @eval begin
    function ($op)(a::Map,b::Map)
      MapFromBroadcastBinaryOp($op,a,b)
    end
  end
end

function varinner(a::Field{D,T},b::Field{D,T}) where {D,T}
  inner(a,b)
end

function varinner(a::Basis{D,T},b::Field{D,T}) where {D,T}
  inner(a,newaxis(b,dim=1))
end

function varinner(a::Basis{D,T},b::Basis{D,T}) where {D,T}
  inner(newaxis(a,dim=2),newaxis(b,dim=1))
end

#"""
#Map generated from a binary operator over two `Map`
#"""
#struct MapFromBinaryOp{O,A,B,S,M,T,N} <: Map{S,M,T,N}
#  op::O
#  a::A
#  b::B
#end
#
#for op in (:+, :-, :*, :/)
#  @eval begin
#    function ($op)(a::Map,b::Map)
#      MapFromBinaryOp($op,a,b)
#    end
#  end
#end
#
#function inner(a::Field{D,T},b::Field{D,T}) where {D,T}
#  MapFromBinaryOp(varinner,a,b)
#end
#
#function inner(a::Basis{D,T},b::Field{D,T}) where {D,T}
#  MapFromBinaryOp(varinner,a,b)
#end
#
#function inner(a::Basis{D,T},b::Basis{D,T}) where {D,T}
#  MapFromBinaryOp(varinner,a,b)
#end
#
#function expand(a::Basis,b::AbstractVector)
#  println("No here")
#  FieldFromExpand(a,b)
#end
#
#
#function MapFromBinaryOp(op::Function,a::Map{S,M,TA,NA},b::Map{S,M,TB,NB}) where {S,M,TA,NA,TB,NB}
#  O = typeof(op)
#  A = typeof(a)
#  B = typeof(b)
#  R = Base._return_type(op,Tuple{TA,TB})
#  N = max(NA,NB)
#  MapFromBinaryOp{O,A,B,S,M,R,N}(op,a,b)
#end
#
#function evaluate(self::MapFromBinaryOp{O,A,B,S,M,T,N},input::AbstractArray{S,M}) where {O,A,B,S,M,T,N}
#  avals = evaluate(self.a,input)
#  bvals = evaluate(self.b,input)
#  self.op.(avals,bvals)
#end
#
#struct FieldFromExpand{
#  D,S,R,T,B<:Basis{D,S},C<:AbstractVector{R}} <: Field{D,T}
#  basis::B
#  coeffs::C
#end
#
#function FieldFromExpand(basis::Basis{D,S},coeffs::AbstractVector{R}) where {D,S,R}
#  T = Base._return_type(outer,Tuple{S,R})
#  B = typeof(basis)
#  C = typeof(coeffs)
#  FieldFromExpand{D,S,R,T,B,C}(basis,coeffs)
#end
#
#function evaluate(self::FieldFromExpand{D},points::AbstractVector{Point{D}}) where D
#  basisvals = evaluate(self.basis,points)
#  expand(basisvals,self.coeffs)
#end
#
#function gradient(self::FieldFromExpand)
#  gradbasis = gradient(self.basis)
#  FieldFromExpand(gradbasis,self.coeffs)
#end
#
#varinner(a::FieldValue,b::FieldValue) = inner(a,b)
#
#const BasisValues{T} = AbstractArray{T,2} where T <: FieldValue
#
#varinner(a::BasisValues{T},b::T) where T<:FieldValue = inner(a,newaxis(b,dim=1))
#
#varinner(a::BasisValues{T},b::BasisValues{T}) where T<:FieldValue = inner(newaxis(a,dim=2),newaxis(b,dim=1))
#
## expand(a::BasisValues,b) = cellsum(outer.(a,newaxis(b,dim=2)),dim=1)
#function expand(a::AbstractArray{T,2}, b::AbstractVector{S}) where {T,S}
#  num_dofs = length(b)
#  @assert size(a)[1] == length(b)
#  num_points = size(a)[end]
#  R = Base._return_type(outer,Tuple{T,S})
#  c = zeros(R,num_points)
#  for i in 1:num_dofs
#    for j in 1:num_points
#      c[j] += a[i,j]*b[i]
#    end
#  end
#  return c
#end
## @santiagobadia : Think about efficiency improvements or other methods to
## implement this function
#
#function newaxis(self::AbstractArray;dim::Int)
#  s = [ v for v in size(self)]
#  insert!(s,dim,1)
#  shape = tuple(s...)
#  return copy(reshape(self,shape))
#end
