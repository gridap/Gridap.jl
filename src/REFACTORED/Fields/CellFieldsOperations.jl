module CellFieldsOperations

using Gridap
using Gridap.Kernels: VarinnerKernel
using Gridap.Kernels: LinCombKernel
using Gridap.Kernels: PhysGradKernel
using Base: @pure
using Base: @propagate_inbounds

export varinner
export lincomb
export attachgeomap
export compose
import Gridap: evaluate
import Gridap: gradient
import Gridap: HasGradientStyle
import Base: iterate
import Base: length
import Base: size
import Base: getindex
import Base: IndexStyle
import Base: +, -
import Base: ∘

for op in (:+,:-)
  @eval begin

    function ($op)(a::CellFieldLike)
      sa = HasGradientStyle(a)
      _compute_sum_or_sub($op,a,sa)
    end

    function ($op)(a::CellFieldLike,b::CellFieldLike)
      sa = HasGradientStyle(a)
      sb = HasGradientStyle(b)
      _compute_sum_or_sub($op,a,b,sa,sb)
    end

  end
end

function varinner(a::CellField,b::CellField)
  inner(a,b)
end

function varinner(a::CellBasis,b::CellFieldLike)
  k = VarinnerKernel()
  apply(k,a,b)
end

function lincomb(a::CellBasis,b::CellVector)
  sa = HasGradientStyle(a)
  _lincomb(a,b,sa)
end

function attachgeomap(a::CellBasis{D},b::CellGeomap{D,D}) where D
  @assert HasGradientStyle(a) == GradientYesStyle()
  refg = gradient(a)
  jac = gradient(b)
  k = PhysGradKernel()
  physg = apply(k,jac,refg)
  _merge_val_and_grad(a,physg)
end

function compose(f::Function,g::CellFieldLike)
  h = hasmethod(gradient,(typeof(f),) )
  _compose(Val(h),f,g)
end

(∘)(f::Function,g::CellFieldLike) = compose(f,g)

function compose(f::Function,w::Vararg{<:CellFieldLike})
  h = hasmethod(gradient,(typeof(f),) )
  _compose(Val(h),f,w...)
end

function _lincomb(a,b,sa::GradientYesStyle)
  k = LinCombKernel()
  v = apply(k,a,b)
  ag = gradient(a)
  g = apply(k,ag,b)
  _merge_val_and_grad(v,g)
end

function _lincomb(a,b,sa)
  k = LinCombKernel()
  apply(k,a,b)
end

function _compose(::Val{true},f,u...)
  fgrad = gradient(f)
  v = apply(f,u...,broadcast=true)
  g = apply(fgrad,u...,broadcast=true)
  _merge_val_and_grad(v,g)
end

function _compose(::Val{false},f,u...)
  v = apply(f,u...,broadcast=true)
  v
end

function _compute_sum_or_sub(op,a,::GradientYesStyle)
  v = apply(op,a,broadcast=true)
  g = apply(op,gradient(a),broadcast=true)
  _merge_val_and_grad(v,g)
end

function _compute_sum_or_sub(op,a,sa)
  v = apply(op,a,broadcast=true)
  v
end

function _compute_sum_or_sub(op,a,b,::GradientYesStyle,::GradientYesStyle)
  v = apply(op,a,b,broadcast=true)
  g = apply(op,gradient(a),gradient(b),broadcast=true)
  _merge_val_and_grad(v,g)
end

function _compute_sum_or_sub(op,a,b,sa,sb)
  v = apply(op,a,b,broadcast=true)
  v
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

HasGradientStyle(::Type{<:IterCellFieldLikeAndGradient}) = GradientYesStyle()

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

HasGradientStyle(::Type{<:IndexCellFieldLikeAndGradient}) = GradientYesStyle()

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

