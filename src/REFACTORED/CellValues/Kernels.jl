module Kernels

using Test
using Gridap
using Gridap.Helpers
using Base.Cartesian: @nloops, @nexprs, @nref
using TensorValues

export NumberKernel
export ArrayKernel
export NumberKernelFromFunction
export ArrayKernelFromBroadcastedFunction
export compute_type
export compute_value
export compute_value!
export compute_size
export compute_ndim
export test_number_kernel
export test_array_kernel

# Interfaces

abstract type NumberKernel end

function compute_type(::NumberKernel,::Vararg{<:Type})::Type{<:NumberLike}
  @abstractmethod
end

function compute_value(::NumberKernel,::Vararg)::NumberLike
  @abstractmethod
end

abstract type ArrayKernel end

function compute_type(::ArrayKernel,::Vararg{<:Type})::Type{<:NumberLike}
  @abstractmethod
end

function compute_ndim(::ArrayKernel,::Vararg{Int})::Int
  @abstractmethod
end

function compute_size(::ArrayKernel,::Vararg{<:NTuple})::NTuple
  @abstractmethod
end

function compute_value!(::AbstractArray,::ArrayKernel,::Vararg)
  @abstractmethod
end

# Testers

function test_number_kernel(k::NumberKernel,o::T,i::Vararg) where T
  t = [ typeof(ii) for ii in i ]
  S = compute_type(k,t...)
  @test S == T
  r = compute_value(k,i...)
  @test r == o
end

function test_array_kernel(
  k::ArrayKernel,o::AbstractArray{T,N},i::Vararg) where {T,N}
  t = [ eltype(ii) for ii in i ]
  S = compute_type(k,t...)
  @test S == T
  s = [ _size_for_broadcast(ii) for ii in i ]
  si = compute_size(k,s...)
  @test size(o) == si
  d = [length(si) for si in s]
  n = compute_ndim(k,d...)
  @test n == length(si)
  r = Array{T,N}(undef,si)
  compute_value!(r,k,i...)
  @test r == o
end

@inline _size_for_broadcast(a) = size(a)

@inline _size_for_broadcast(a::NumberLike) = ()

# Implementations

struct NumberKernelFromFunction{F<:Function} <: NumberKernel
  fun::F
end

function compute_type(
  self::NumberKernelFromFunction,types::Vararg{<:Type})
  T = Base._return_type(self.fun,types)
  @assert T <: NumberLike
  T
end

function compute_value(self::NumberKernelFromFunction,args::Vararg)
  self.fun(args...)
end

struct ArrayKernelFromBroadcastedFunction{F<:Function} <: ArrayKernel
  fun::F
end

function compute_type(
  self::ArrayKernelFromBroadcastedFunction,etypes::Vararg{<:Type})
  Base._return_type(self.fun,etypes)
end

function compute_size(::ArrayKernelFromBroadcastedFunction,s::Vararg{<:NTuple})
  Base.Broadcast.broadcast_shape(s...)
end

function compute_value!(
  v::AbstractArray, k::ArrayKernelFromBroadcastedFunction, a::Vararg)
  broadcast!(k.fun,v,a...)
end

function compute_ndim(::ArrayKernelFromBroadcastedFunction,d::Vararg{Int})
  n = 0
  for di in d
    n = max(n,di)
  end
  n
end

struct CellSumKernel{D} <: ArrayKernel end

compute_type(k::CellSumKernel,T::Type) = T

compute_ndim(k::CellSumKernel,nd::Int) = nd-1

@generated function compute_size(
  k::CellSumKernel{D},asize::NTuple{N,Int}) where {D,N}

  @assert N > 0
  @assert D <= N
  str = join([ "asize[$i]," for i in 1:N if i !=D ])
  Meta.parse("($str)")
end

@generated function compute_value!(
  B::AbstractArray,
  ::CellSumKernel{D},
  A::AbstractArray{T,N}) where {T,N,D}

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

struct LinCombKernel <: ArrayKernel end

function compute_type(::LinCombKernel,::Type{T},::Type{S}) where {T,S}
  Base._return_type(outer,Tuple{T,S})
end

function compute_ndim(::LinCombKernel,n::Int,m::Int)
  @assert n == 2
  @assert m == 1
  1
end

function compute_size(::LinCombKernel,sa::NTuple{2,Int},sb::NTuple{1,Int})
  ndofs, npoints = sa
  ndofsb, = sb
  @assert ndofsb == ndofs
  (npoints,)
end

function compute_value!(
  v::AbstractArray{T},::LinCombKernel,a::AbstractArray,b::AbstractArray) where T
  ndofs, npoints = size(a)
  for i in eachindex(v)
    @inbounds v[i] = zero(T)
  end
  for j in 1:npoints
    for i in 1:ndofs
      @inbounds v[j] += outer(a[i,j],b[i])
    end
  end
end

struct VarinnerKernel <: ArrayKernel end

function compute_type(::VarinnerKernel,::Type{T},::Type{S}) where {T,S}
  Base._return_type(inner,Tuple{T,S})
end

function compute_ndim(::VarinnerKernel,n::Int,m::Int)
  (n-1) + (m-1) + 1
end

function compute_size(::VarinnerKernel,sa::NTuple{2,Int},sb::NTuple{1,Int})
  ndofs, npoints = sa
  npointsb, = sb
  @assert npoints == npointsb
  sa
end

function compute_size(::VarinnerKernel,sa::NTuple{2,Int},sb::NTuple{2,Int})
  ndofs, npoints = sa
  ndofsb, npointsb = sb
  @assert npoints == npointsb
  (ndofs,ndofsb,npoints)
end

function compute_value!(
  v::AbstractArray,::VarinnerKernel,a::AbstractArray{T,2},b::AbstractArray{S,1}) where {T,S}
  ndofs, npoints = size(a)
  for j in 1:npoints
    for i in 1:ndofs
      @inbounds v[i,j] = inner(a[i,j],b[j])
    end
  end
end

function compute_value!(
  v::AbstractArray,::VarinnerKernel,a::AbstractArray{T,2},b::AbstractArray{S,2}) where {T,S}
  ndofsa, npoints = size(a)
  ndofsb, _ = size(b)
  for k in 1:npoints
    for j in 1:ndofsb
      for i in 1:ndofsa
        @inbounds v[i,j,k] = inner(a[i,k],b[j,k])
      end
    end
  end
end

end # module
