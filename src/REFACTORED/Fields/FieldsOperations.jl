module FieldsOperations

using Gridap

import Base: +, -

struct FieldLikeAndGradient{D,T,N,X,V,G} <: FieldLike{D,T,N,X}
  field::V
  grad::G
end

function FieldLikeAndGradient(
  field::FieldLike{D,T,N,X},
  grad::FieldLike{D,TG,N,X}) where {D,T,TG,N,X}
  V = typeof(field)
  G = typeof(grad)
  FieldLikeAndGradient{D,T,N,X,V,G}(field,grad)
end

function evaluate!(
  this::FieldLikeAndGradient{D,T,N,X},
  points::AbstractVector{Point{D,X}},
  v::AbstractArray{T,N}) where {D,X,T,N}
  evaluate!(this.field,points,v)
end

function return_size(this::FieldLikeAndGradient,s::Tuple{Int})
  return_size(this.Field,s)
end

gradient(f::FieldLikeAndGradient) = f.grad

end # module
