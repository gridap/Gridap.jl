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
      IterCellFieldLikeAndGradient(v,g)
    end

    function ($op)(f::IndexCellFieldLike)
      v = apply($op,f,broadcast=true)
      g = apply($op,gradient(f),broadcast=true)
      IndexCellFieldLikeAndGradient(v,g)
    end

    function ($op)(a::CellFieldLike,b::CellFieldLike)
      v = apply($op,a,b,broadcast=true)
      g = apply($op,gradient(a),gradient(b),broadcast=true)
      IterCellFieldLikeAndGradient(v,g)
    end

    function ($op)(a::IndexCellFieldLike,b::IndexCellFieldLike)
      v = apply($op,a,b,broadcast=true)
      g = apply($op,gradient(a),gradient(b),broadcast=true)
      IndexCellFieldLikeAndGradient(v,g)
    end

  end
end

struct IterCellFieldLikeAndGradient{D,T,N,R,V,G} <: IterCellFieldLike{D,T,N,R}
  val::V
  grad::G
end

function IterCellFieldLikeAndGradient(
  val::IterCellFieldLike{D,T,N,R},
  grad::CellFieldLike{D,TG,N}) where {D,T,TG,N,R}
  @show T
  @assert length(val) == length(grad)
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

struct IndexCellFieldLikeAndGradient{D,T,N,C,R,V,G} <: IndexCellFieldLike{D,T,N,C,R}
  val::V
  grad::G
end

function IndexCellFieldLikeAndGradient(
  val::IndexCellFieldLike{D,T,N,C,R},
  grad::IndexCellFieldLike{D,TG,N}) where {D,T,TG,N,C,R}
  @assert length(val) == length(grad)
  V = typeof(val)
  G = typeof(grad)
  IndexCellFieldLikeAndGradient{D,T,N,C,R,V,G}(val,grad)
end

gradient(f::IndexCellFieldLikeAndGradient) = f.grad

function evaluate(cm::IndexCellFieldLikeAndGradient{D},ca::CellPoints{D}) where D
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

