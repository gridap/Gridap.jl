module Maps

using Test
using Gridap
using Gridap.Helpers

export Map
export evaluate!
export evaluate
export return_size
export test_map

"""
Abstract map that takes an `M`-dim array of `S` values and returns an `N`-dim
array of `T` values.
For efficiency reasons, `T` should be a concrete type. However, `S` can be also an abstract type or 
a union of types
"""
abstract type Map{S,M,T,N} end

"""
Evaluate a `Map` on a set of points
"""
function evaluate!(
  this::Map{S,M,T,N},
  points::AbstractArray{<:S,M},
  v::AbstractArray{T,N}) where {S,M,T,N}
  @abstractmethod
end

"""
Return dimension of the output array
"""
function return_size(
  ::Map{S,M,T,N},::NTuple{M,Int})::NTuple{N,Int} where {S,M,T,N}
  @abstractmethod
end

"""
Same as `evaluate!` but allocates output
"""
function evaluate(
  this::Map{S,M,T,N},
  points::AbstractArray{<:S,M}) where {S,M,T,N}
  v_size = return_size(this, size(points))
  v = Array{T,N}(undef, v_size)
  evaluate!(this,points,v)
  return v
end

# Testers

function test_map(
  m::Map{S,M,T,N},x::AbstractArray{<:S,M},y::AbstractArray{T,N}) where {S,M,T,N}

  z = evaluate(m,x)

  @test z ≈ y

  s = return_size(m,size(x))

  @test s == size(y)

end

end # module Maps
