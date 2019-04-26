module Maps

using Numa.Helpers
using Numa.FieldValues

export Map
export Field
export AnalyticalField
export GeoMap
export Basis

import Numa: evaluate, gradient
export evaluate!
export return_size

import Base: +, -, *, /, ∘

import Numa.FieldValues: inner, outer

"""
Abstract map that takes an `M`-dim array of `S` values and returns an `N`-dim
array of `T` values. For vectorization purposes, a function that takes `S` and
returns `T` is declared as `Map{S,1,T,1}`
"""
abstract type Map{S,M,T,N} end

"""
Abstract field of rank `T` (e.g., scalar, vector, tensor) on a manifold of
dimension `D`
"""
const Field{D,T} = Map{Point{D},1,T,1} where {D,T<:FieldValue}

"""
Abstract basis for a space of fields of rank `T` (e.g., scalar, vector, tensor)
on a manifold of dimension `D`
"""
const Basis{D,T} = Map{Point{D},1,T,2} where {D,T<:FieldValue}

"""
Abstract geometry map
"""
const Geomap{D,Z} = Field{D,Point{Z}}

"""
Evaluate a `Map` on a set of points
"""
function evaluate!(
  this::Map{S,M,T,N},
  points::AbstractArray{S,M},
  v::AbstractArray{T,N}) where {S,M,T,N}
  @abstractmethod
end

"""
Creates the gradient of a `Map`
"""
function gradient(
  this::Map{S,M,T,N})::Map{S,M,TG,N} where {S,M,T,N,TG}
  @abstractmethod
end

"""
Same as `evaluate!` but allocates output
"""
function evaluate(
  this::Map{S,M,T,N},
  points::AbstractArray{S,M}) where {S,M,T,N}
  v_size = return_size(this, size(points))
  v = Array{T,N}(undef, v_size)
  evaluate!(this,points,v)
  return v
end

"""
Return dimension of the output array
"""
return_size(::Map,::NTuple{N,Int} where N)::Tuple = @abstractmethod

# Concrete structs

"""
Field generated from an analytical function
"""
struct AnalyticalField{D,T,F<:Function} <: Field{D,T}
  fun::F
end

function AnalyticalField(fun::Function,D::Int)
  T = Base._return_type(fun,Tuple{Point{D}})
  F = typeof(fun)
  AnalyticalField{D,T,F}(fun)
end

function evaluate!(
  this::AnalyticalField{D,T},
  points::AbstractVector{Point{D}},
  v::AbstractVector{T}) where {D,T}
  v .= this.fun.(points)
end

function gradient(this::AnalyticalField{D}) where D
  gradfun = gradient(this.fun)
  AnalyticalField(gradfun,D)
end
return_size(::AnalyticalField, p_size::Tuple{Vararg{Int64,N}} where N) = p_size
return_size(::AnalyticalField, p_size::Int) = (p_size)

include("Operators.jl")

end # module Maps
