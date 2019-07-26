module CellFieldsOperations

using Gridap
using Gridap.CachedValues
using Gridap.Kernels: VarinnerKernel
using Gridap.Kernels: LinCombKernel
using Gridap.Kernels: PhysGradKernel
using Gridap.FieldsOperations: FieldLikeAndGradient
using Base: @pure
using Base: @propagate_inbounds

using Gridap.FieldsOperations: _gradient
using Gridap.CellMapApply: CellMapFromKernel, IndexCellMapFromKernel

export varinner
export lincomb
export attachgeomap
export compose
export symmetric_gradient
export ε
import Gridap: evaluate
import Gridap: gradient
import Gridap: reindex
import Base: iterate
import Base: length
import Base: size
import Base: getindex
import Base: IndexStyle
import Base: ∘

function varinner(a::CellField,b::CellField)
  inner(a,b)
end

function varinner(a::CellBasis,b::CellFieldLike)
  k = VarinnerKernel()
  apply(k,a,b)
end

function lincomb(a::CellBasis,b::CellVector)
  k = LinCombKernel()
  v = apply(k,a,b)
  ag = gradient(a)
  g = apply(k,ag,b)
  _merge_val_and_grad(v,g)
end

function attachgeomap(a::CellBasis{D},b::CellGeomap{D,D}) where D
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

function symmetric_gradient(f::CellFieldLike{D,T}) where {D,T<:VectorValue}
  g = gradient(f)
  apply(symmetic_part,g,broadcast=true)
end

const ε = symmetric_gradient

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
  R = FieldLikeAndGradient{D,T,N,eltype(val),eltype(grad)}
  V = typeof(val)
  G = typeof(grad)
  IterCellFieldLikeAndGradient{D,T,N,R,V,G}(val,grad)
end

gradient(f::IterCellFieldLikeAndGradient) = f.grad

function evaluate(cm::IterCellFieldLikeAndGradient{D},ca::CellPoints{D}) where D
  evaluate(cm.val,ca)
end

@inline function iterate(f::IterCellFieldLikeAndGradient)
  fnext = iterate(f.val)
  gnext = iterate(f.grad)
  if fnext === nothing; return nothing; end
  if gnext === nothing; return nothing; end
  fi, fstate = fnext
  gi, gstate = gnext
  v = FieldLikeAndGradient(fi,gi)
  state = (v,fstate,gstate)
  (v,state)
end

@inline function iterate(f::IterCellFieldLikeAndGradient,state)
  v, fstate, gstate = state
  fnext = iterate(f.val,fstate)
  gnext = iterate(f.grad,gstate)
  if fnext === nothing; return nothing; end
  if gnext === nothing; return nothing; end
  fi, fstate = fnext
  gi, gstate = gnext
  v.val = fi
  v.grad = gi
  state = (v,fstate,gstate)
  (v,state)
end

length(f::IterCellFieldLikeAndGradient) = length(f.val)

struct IndexCellFieldLikeAndGradient{
  D,T,N,C,R<:FieldLike{D,T,N},V,G,F} <: IndexCellValue{R,C}
  val::V
  grad::G
  cache::CachedValue{F}
end

function IndexCellFieldLikeAndGradient(
  val::IndexCellFieldLike{D,T,N,C},
  grad::IndexCellFieldLike{D,TG,N}) where {D,T,TG,N,C}
  @assert length(val) == length(grad)
  R = FieldLikeAndGradient{D,T,N,eltype(val),eltype(grad)}
  F = Union{Nothing,R}
  V = typeof(val)
  G = typeof(grad)
  cache = CachedValue{F}(nothing)
  IndexCellFieldLikeAndGradient{D,T,N,C,R,V,G,F}(val,grad,cache)
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
  fi = f.val[i...]
  gi = f.grad[i...]
  if f.cache.value === nothing
    f.cache.value = FieldLikeAndGradient(fi,gi)
  else
    f.cache.value.val = fi
    f.cache.value.grad = gi
  end
  f.cache.value
end

function gradient(cm::CellMapFromKernel)
  _gradient(cm.kernel,cm.cellvalues...)
end

function gradient(cm::IndexCellMapFromKernel)
  _gradient(cm.kernel,cm.cellvalues...)
end

function reindex(cm::IterCellFieldLikeAndGradient,indices::CellValue{<:IndexLike})
  _reindex(cm,indices)
end

function reindex(cm::IndexCellFieldLikeAndGradient,indices::CellValue{<:IndexLike})
  _reindex(cm,indices)
end

function reindex(cm::IndexCellFieldLikeAndGradient,indices::IndexCellValue{<:IndexLike})
  _reindex(cm,indices)
end

function _reindex(cm,indices)
  v = reindex(cm.val,indices)
  g = reindex(cm.grad,indices)
  _merge_val_and_grad(v,g)
end

end # module
