module FieldsOperations

using Gridap
using Gridap.Helpers

export change_basis
import Gridap: evaluate!
import Gridap: return_size
import Gridap: gradient
import Gridap.MapApply: MapFromKernel

function change_basis(basis::Basis,changeofbasis::Matrix)
  BasisFromChangeOfBasis(basis,changeofbasis)
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

struct BasisFromChangeOfBasisGrad{D,T,C,B<:Basis} <: Basis{D,T}
  basis::B
  changeofbasis::Matrix{C}
  v::Vector{T}
end

function BasisFromChangeOfBasisGrad(
  basis::Basis{D,T},changeofbasis::Matrix{C}) where {D,T,C}
  B = typeof(basis)
  S = Base._return_type(*,(C,T))
  v = zeros(S,length(basis))
  BasisFromChangeOfBasisGrad{D,S,C,B}(basis,changeofbasis,v)
end

struct BasisFromChangeOfBasis{D,T,C,B<:Basis,G} <: Basis{D,T}
  basis::B
  changeofbasis::Matrix{C}
  g::G
  v::Vector{T}
end

function BasisFromChangeOfBasis(
  basis::Basis{D,T},changeofbasis::Matrix{C}) where {D,T,C}
  B = typeof(basis)
  g = BasisFromChangeOfBasisGrad(gradient(basis),changeofbasis)
  G = typeof(g)
  S = Base._return_type(*,(C,T))
  v = zeros(S,length(basis))
  BasisFromChangeOfBasis{D,S,C,B,G}(basis,changeofbasis,g,v)
end

const _BasisFromChangeOfBasis{D,T} = Union{
  BasisFromChangeOfBasis{D,T}, BasisFromChangeOfBasisGrad{D,T}}

function evaluate!(
  this::_BasisFromChangeOfBasis{D,T},
  points::AbstractVector{<:Point{D}},
  v::AbstractArray{T,N}) where {D,T,N}
  evaluate!(this.basis,points,v)
  ndofs, npoints = size(v)
  for j in 1:npoints
    for i in 1:ndofs
      a = zero(T)
      for k in 1:ndofs
        @inbounds a += this.changeofbasis[i,k]*v[k,j]
      end
      this.v[i] = a
    end
    for i in 1:ndofs
      @inbounds v[i,j] = this.v[i]
    end
  end
end

function return_size(
  f::_BasisFromChangeOfBasis,s::NTuple{N,Int} where N)
  return_size(f.basis,s)
end

gradient(f::BasisFromChangeOfBasis) = f.g

gradient(f::BasisFromChangeOfBasisGrad) = @notimplemented

function gradient(m::MapFromKernel)
  _gradient(m.kernel,m.inputs...)
end

_gradient(k,ms...) = @notimplemented

for op in (:+, :-)
  @eval begin
    function _gradient(k::ArrayKernelFromBroadcastedFunction{typeof($op)},ms...)
      g = [ gradient(m) for m in ms  ]
      apply(k,g...)
    end
  end
end

end # module
