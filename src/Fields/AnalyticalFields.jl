module AnalyticalFields

using Gridap
export AnalyticalField

import Gridap: evaluate!
import Gridap: return_size
import Gridap.FieldsOperations: FieldLikeAndGradient

"""
Field generated from an analytical function
"""
struct AnalyticalField{D,T,F<:Function} <: Field{D,T}
  fun::F
end

function AnalyticalField(fun::Function,D::Int,X::Type=Float64)
  h = hasmethod(gradient,(typeof(fun),) )
  _setup(Val(h),fun,D,X)
end

function _setup(::Val{true},fun,D,X)
  gun = gradient(fun)
  v = _AnalyticalField(fun,D,X)
  g = _AnalyticalField(gun,D,X)
  FieldLikeAndGradient(v,g)
end

function _setup(::Val{false},fun,D,X)
  v = _AnalyticalField(fun,D,X)
  v
end

function _AnalyticalField(fun::Function,D::Int,X::Type)
  F = typeof(fun)
  T = Base._return_type(fun,Tuple{Point{D,X}})
  AnalyticalField{D,T,F}(fun)
end

function evaluate!(
  this::AnalyticalField{D,T},
  points::AbstractVector{<:Point{D}},
  v::AbstractVector{T}) where {D,T}
  v .= this.fun.(points)
end

return_size(::AnalyticalField, psize::NTuple{N,Int} where N) = psize

end # module
