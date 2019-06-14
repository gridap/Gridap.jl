module FieldsOperations

using Gridap

import Gridap: evaluate!
import Gridap: return_size
import Gridap: gradient
import Base: +, -

for op in (:+,:-)
  @eval begin

    function ($op)(f::FieldLike)
      v = apply($op,f,broadcast=true)
      g = apply($op,gradient(f),broadcast=true)
      FieldLikeAndGradient(v,g)
    end

    function ($op)(a::FieldLike,b::FieldLike)
      v = apply($op,a,b,broadcast=true)
      g = apply($op,gradient(a),gradient(b),broadcast=true)
      FieldLikeAndGradient(v,g)
    end

  end
end

mutable struct FieldLikeAndGradient{D,T,N,V,G} <: FieldLike{D,T,N}
  val::V
  grad::G
end

function FieldLikeAndGradient(
  val::FieldLike{D,T,N},
  grad::FieldLike{D,TG,N}) where {D,T,TG,N}
  V = typeof(val)
  G = typeof(grad)
  FieldLikeAndGradient{D,T,N,V,G}(val,grad)
end

function evaluate!(
  this::FieldLikeAndGradient{D,T,N},
  points::AbstractVector{<:Point{D}},
  v::AbstractArray{T,N}) where {D,T,N}
  evaluate!(this.val,points,v)
end

function return_size(
  f::FieldLikeAndGradient,s::NTuple{N,Int} where N)
  return_size(f.val,s)
end

gradient(f::FieldLikeAndGradient) = f.grad

end # module
