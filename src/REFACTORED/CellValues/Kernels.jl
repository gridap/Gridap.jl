module Kernels

using Test
using Gridap
using Gridap.Helpers

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

function compute_value!(::AbstractArray,::NumberKernel,::Vararg)
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

end # module Kernels
