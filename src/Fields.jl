module Fields

using Numa.Helpers
using Numa.FieldValues

export Field
export AnalyticalField

export evaluate
export evaluate!
export gradient

"""
Abstract field of rank `T` (e.g., scalar, vector, tensor) on a manifold of
dimension `D`
"""
abstract type Field{D,T<:FieldValue} end

"""
Evaluate the field on a set of points
"""
function evaluate!(
  this::Field{D,T},
  points::AbstractVector{Point{D}},
  v::AbstractVector{T}) where {D,T<:FieldValue}
  @abstractmethod
end

"""
Creates the gradient field of a field
"""
function gradient(
  this::Field{D,T})::Field{D,TG} where {D,T<:FieldValue,TG<:FieldValue}
  @abstractmethod
end

"""
Same as evaluate! but allocates output
"""
function evaluate(
  this::Field{D,T},points::AbstractVector{Point{D}}) where {D,T<:FieldValue}
  v = Array{T,1}(undef, (length(points),) )
  evaluate!(this,points,v)
  v
end

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
  v::AbstractVector{T}) where {D,T<:FieldValue}
  v .= this.fun.(points)
end

function gradient(this::AnalyticalField{D}) where D
  gradfun = gradient(this.fun)
  AnalyticalField(gradfun,D)
end

end # module Fields
