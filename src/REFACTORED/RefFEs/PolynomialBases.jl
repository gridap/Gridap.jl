module PolynomialBases

using Gridap
using Gridap.Helpers
import TensorPolynomialBases; const tp = TensorPolynomialBases

export MonomialBasis
export PolynomialBasis
import Gridap: evaluate!
import Gridap: return_size
import Gridap: gradient

struct PolynomialBasisGrad{D,T,B,C} <: Basis{D,T}
  basis::B
  cache::C
  v::Vector{T}
end

struct PolynomialBasis{D,T,TG,B,C} <: Basis{D,T}
  basis::B
  cache::C
  v::Vector{T}
  g::PolynomialBasisGrad{D,TG,B,C}
end

"""
Construct a PolynomialBasis formed by monomials up to a given order.
Fine control of which monomials form the basis is done with the filter function.
"""
function MonomialBasis(D::Int,T::Type,filter::Function,order::Int)
  P = Point{D}
  basis = tp.MonomialBasis{P,T}(filter,order)
  PolynomialBasis(basis)
end

"""
Construct a PolynomialBasis formed my the monomials of an anisotropic Q-space
"""
function MonomialBasis(T::Type,orders::NTuple{D,Int}) where D
  P = Point{D}
  basis = tp.MonomialBasis{P,T}(orders)
  PolynomialBasis(basis)
end

function MonomialBasis(T::Type,orders::Vector{Int})
  MonomialBasis(T,tuple(orders...))
end

"""
Construct a PolynomialBasis from a given TensorPolynomialBasis of the
TensorPolynomialBases package.
This allows one to use any of the bases defined in the TensorPolynomialBases
package within the Gridap Project.
"""
function PolynomialBasis(basis::tp.TensorPolynomialBasis)
  cache = tp.ScratchData(basis)
  T = tp.value_type(basis)
  TG = tp.gradient_type(basis)
  D = ndims(basis)
  ndofs = length(basis)
  B = typeof(basis)
  C = typeof(cache)
  v = zeros(T,ndofs)
  vg = zeros(TG,ndofs)
  g = PolynomialBasisGrad{D,TG,B,C}(basis,cache,vg)
  PolynomialBasis{D,T,TG,B,C}(basis,cache,v,g)
end

function evaluate!(
  this::PolynomialBasis{D,T},
  points::AbstractVector{<:Point{D}},
  v::AbstractMatrix{T}) where {D,T}
  ndofs = length(this.v)
  for j in eachindex(points)
    p = points[j]
    tp.evaluate!(this.v,this.basis,p,this.cache)
    for i in 1:ndofs
      v[i,j] = this.v[i]
    end
  end
end

function evaluate!(
  this::PolynomialBasisGrad{D,T},
  points::AbstractVector{<:Point{D}},
  v::AbstractMatrix{T}) where {D,T}
  ndofs = length(this.v)
  for j in eachindex(points)
    p = points[j]
    tp.gradient!(this.v,this.basis,p,this.cache)
    for i in 1:ndofs
      v[i,j] = this.v[i]
    end
  end
end

function return_size(
  this::Union{PolynomialBasis,PolynomialBasisGrad},
  s::NTuple{N,Int} where N)
  ndofs = length(this.basis)
  npoints, = s
  (ndofs,npoints)
end

gradient(f::PolynomialBasis) = f.g

gradient(f::PolynomialBasisGrad) = @notimplemented

end # module
