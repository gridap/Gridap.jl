module Maps

using Numa.Helpers
using Numa.FieldValues

export Map
export Field
export AnalyticalField

import Numa: evaluate, gradient
export evaluate!

"""Create a new object Map{S,M,T,N} that given an M-dim array of S provides an N-dim array of T
Abstract map that takes an `M`-dim array of `S` values and returns an `N`-dim array of `T` values
"""
abstract type Map{S,M,T,N} end

"""
Abstract field of rank `T` (e.g., scalar, vector, tensor) on a manifold of
dimension `D`
"""
const Field{D,T} = Map{Point{D},1,T,1}

"""
Abstract basis for a field of rank `T` (e.g., scalar, vector, tensor) on a manifold of
dimension `D`
"""
# cont Basis{D,T} = Map{Point{D},1,T,2}

"""
Evaluate the Map on a set of points
"""
function evaluate!(
  this::Map{S,M,T,N},
  points::AbstractVector{AbstractArray{S,M}},
  v::AbstractVector{Array{T,N}}) where {S,M,T,N}
  @abstractmethod
end

"""
Creates the gradient Map of a Map
"""
function gradient(
  this::Map{S,M,T,N})::Map{S,M,TG,N} where {S,M,T,N,TG}
  @abstractmethod
end

"""
Same as evaluate! but allocates output
"""
function evaluate(
  this::Map{S,M,T,N},
  points::AbstractVector{AbstractArray{S,M}}) where {S,M,T,N}
  @abstractmethod
end

"""
evaluate! for `Field`
"""
function evaluate(
  this::Field{D,T},points::AbstractVector{Point{D}}) where {D,T<:FieldValue}
  v = Array{T,1}(undef, (length(points),) )
  evaluate!(this,points,v)
  v
end
# @santiagobadia: Not sure this is the way to go. Should it be in Map, it does
# not seem posible, since not concrete the way AbstractArray{T,N} stored,
# abstract at this point

"""
Field generated from an analytical function
"""
struct AnalyticalField{D,T<:FieldValue,F<:Function} <: Field{D,T}
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

end # module Maps
