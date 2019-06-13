module CellFieldsOperations

using Gridap
using Base: @pure
using Base: @propagate_inbounds

import Gridap: evaluate
import Gridap: gradient
import Base: iterate
import Base: length
import Base: size
import Base: getindex
import Base: IndexStyle
import Base: +, -

for op in (:+,:-)
  @eval begin

    function ($op)(f::CellFieldLike)
      v = apply($op,f,broadcast=true)
      g = apply($op,gradient(f),broadcast=true)
      _merge_val_and_grad(v,g)
    end

    function ($op)(a::CellFieldLike,b::CellFieldLike)
      v = apply($op,a,b,broadcast=true)
      g = apply($op,gradient(a),gradient(b),broadcast=true)
      _merge_val_and_grad(v,g)
    end

  end
end

function _merge_val_and_grad(v::CellFieldLike,g::CellFieldLike)
  IterCellFieldLikeAndGradient(v,g)
end

function _merge_val_and_grad(v::IndexCellFieldLike,g::IndexCellFieldLike)
  IndexCellFieldLikeAndGradient(v,g)
end

struct IterCellFieldLikeAndGradient{D,T,N,R<:FieldLike{D,T,N},V,G} <: IterCellValue{R}
  val::V
  grad::G
end

function IterCellFieldLikeAndGradient(
  val::IterCellFieldLike{D,T,N}, grad::CellFieldLike{D,TG,N}) where {D,T,TG,N}
  @assert length(val) == length(grad)
  R = eltype(val)
  V = typeof(val)
  G = typeof(grad)
  IterCellFieldLikeAndGradient{D,T,N,R,V,G}(val,grad)
end

gradient(f::IterCellFieldLikeAndGradient) = f.grad

function evaluate(cm::IterCellFieldLikeAndGradient{D},ca::CellPoints{D}) where D
  evaluate(cm.val,ca)
end

@inline iterate(f::IterCellFieldLikeAndGradient) = iterate(f.val)

@inline iterate(f::IterCellFieldLikeAndGradient,state) = iterate(f.val,state)

length(f::IterCellFieldLikeAndGradient) = length(f.val)

struct IndexCellFieldLikeAndGradient{
  D,T,N,C,R<:FieldLike{D,T,N},V,G} <: IndexCellValue{R,C}
  val::V
  grad::G
end

function IndexCellFieldLikeAndGradient(
  val::IndexCellFieldLike{D,T,N,C},
  grad::IndexCellFieldLike{D,TG,N}) where {D,T,TG,N,C}
  @assert length(val) == length(grad)
  R = eltype(val)
  V = typeof(val)
  G = typeof(grad)
  IndexCellFieldLikeAndGradient{D,T,N,C,R,V,G}(val,grad)
end

gradient(f::IndexCellFieldLikeAndGradient) = f.grad

function evaluate(
  cm::IndexCellFieldLikeAndGradient{D},ca::CellPoints{D}) where D
  evaluate(cm.val,ca)
end

size(f::IndexCellFieldLikeAndGradient) = size(f.val)

@pure function IndexStyle(
  ::Type{IndexCellFieldLikeAndGradient{D,T,N,C,R,V,G}}) where {D,T,N,C,R,V,G}
  IndexStyle(V)
end

@propagate_inbounds function getindex(
  f::IndexCellFieldLikeAndGradient,i::Vararg{<:Integer})
  f.val[i...]
end

end # module

