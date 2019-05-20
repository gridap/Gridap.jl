module Polynomials

using Gridap #@fverdugo to be eliminated
using Gridap.Helpers
using Gridap.FieldValues
using Gridap.Maps: Basis
using Gridap.FieldValues: mutable
using StaticArrays: MVector
# using Gridap.Quadratures
# using Base.Cartesian

export PolynomialBasis
export TensorProductMonomialBasis
export PolynomialBasisWithChangeOfBasis
# export UnivariatePolynomialBasis
# export UnivariateMonomialBasis

import Gridap: evaluate
import Gridap: gradient, ∇
import Gridap.Maps: evaluate!
import Gridap.Maps: return_size

"""
Abstract type representing a multivariate polynomial basis
with value of type T in a coordinate space of D dimensions
"""
abstract type PolynomialBasis{D,T} <: Basis{D,T} end

Base.length(::PolynomialBasis{D,T} where {D,T})::Int = @abstractmethod
return_size(this::PolynomialBasis, p_size::NTuple{N,Int} where N) = (length(this),p_size...)

"""
Same as evaluate! but allocates output
"""
function evaluate(self::PolynomialBasis{D,T},points::AbstractArray{Point{D},1}) where {D,T}
  vals = Array{T,2}(undef,(length(self),length(points)))
  evaluate!(self,points,vals)
  vals
end

"""
First axis of v for dofs, second for points
"""
evaluate!(::PolynomialBasis{D,T},::AbstractVector{Point{D}},v::AbstractArray{T,2}) where {D,T} = @abstractmethod

"""
Object that represents the gradient of the elements of a basis
"""
struct GradPolynomialBasis{D,TG,T,B<:PolynomialBasis{D,T}} <: PolynomialBasis{D,TG}
  basis::B
end

"""
evaluate! overwritten for gradients of bases
"""
function evaluate!(self::GradPolynomialBasis{D,T},
  points::AbstractArray{Point{D},1}, v::AbstractArray{T,2}) where {D,T}
  evaluategradients!(self.basis,points,v)
end

"""
Create a `GradPolynomialBasis` object from a basis
`PolynomialBasis`. The result is a
PolynomialBasis{TG,D} where TG is a type whose rank is one unit
greater than the one of T
"""
function gradient(self::PolynomialBasis{D,T}) where{D,T}
  TG = outer(Point{D},T)
  # TG = Base._return_type(outer,Tuple{Point{D},T})
  B = typeof(self)
  GradPolynomialBasis{D,TG,T,B}(self)
end

Base.length(this::GradPolynomialBasis)::Int = length(this.basis)

function gradient(::GradPolynomialBasis{D,T,B} where{D,T,B})
 error("Gradient of a gradient not available")
end


"""
Abstract type representing a univariate polynomial basis in dimension one
"""
abstract type UnivariatePolynomialBasis <: PolynomialBasis{1,ScalarValue} end

"""
Univariate monomial basis of a given `order`
"""
struct UnivariateMonomialBasis <: UnivariatePolynomialBasis
  order::Int64
end

Base.length(this::UnivariateMonomialBasis)::Int = this.order+1

"""
Auxiliary function that does the same as evaluate but using a pre-allocated
array
"""
function evaluate!(this::UnivariateMonomialBasis,
  points::AbstractVector{Point{1}},v::AbstractArray{ScalarValue,2})
  for (j,p) in enumerate(points)
    for i in 1:length(this)
      v[i,j] = p[1]^(i-1)
    end
  end
end

function evaluategradients!(this::UnivariateMonomialBasis,
  points::AbstractVector{Point{1}},v::AbstractArray{VectorValue{1},2})
  # numd = 1 # Changing this number, we get higher order derivatives, to be used
  # for Hessians, etc.
  for (j,p) ∈ enumerate(points)
    v[1,j] = VectorValue{1}(0.0)
    for i in 2:length(this)
      # val = (i<=numd) ? 0.0 : prod([i-k-1 for k=0:numd-1])p[1]^(i-numd-1)
      v[i,j] = VectorValue{1}((i-1)*p[1]^(i-2))
    end
  end
end

function evaluategradients(self::UnivariateMonomialBasis,
  points::AbstractVector{Point{1}})
  vals = Array{VectorValue{1},2}(undef,(length(self),length(points)))
  evaluategradients!(self,points,vals)
  vals
end

"""
Multivariate monomial basis obtained as tensor product of univariate polynomial
basis per dimension
"""
struct TensorProductMonomialBasis{D,T} <: PolynomialBasis{D,T}
  univariatebases::NTuple{D,UnivariateMonomialBasis}
end

"""
Provide a `TensorProductMonomialBasis` for a vector `order` providing the order per
dimension
"""
function TensorProductMonomialBasis{D,T}(order::Vector{Int64}) where {D,T}
  @assert(length(order) == D)
  uvmbs = [UnivariateMonomialBasis(order[i]) for i=1:length(order)]
  TensorProductMonomialBasis{D,T}(tuple(uvmbs...))
end

function Base.length(this::TensorProductMonomialBasis{D,T})::Int where {D,T}
  p = 1
  for i in 1:D
    p *= length(this.univariatebases[i])
  end
  length(T)*p
end

function evaluate!(this::TensorProductMonomialBasis{D,T},
  points::AbstractVector{Point{D}}, v::AbstractArray{T,2}) where {D,T}
  tpcoor = i -> [ Point{1}(p[i]) for p in points]
  cooruv = [tpcoor(i) for i in 1:D]
  univals = [evaluate(this.univariatebases[i],cooruv[i]) for i in 1:D]
  # @santiagobadia : In the future, we can create a new evaluate! interface with
  # an additional scratch data with univals. To be used in unfitted FEM
  # @santiagobadia : Strided array of points foer every dim instead of cooruv
  cid = ntuple(i -> 1:length(this.univariatebases[i]), D)
  lent = length(T)
  cid = (cid..., 1:lent)
  cid = CartesianIndices(cid)
  E = eltype(T)
  MT = mutable(T)
  aux = zero(MT)
  for (i,j) in enumerate(cid)
    d = j[D+1]
    for k in 1:length(points)
      val = 1.0
      for l in 1:D
        val *= univals[l][j[l],k]
      end
      v[i,k] = insertentry!(val,aux,d)
      insertentry!(zero(E),aux,d)
    end
  end
end

insertentry!(a::Float64, b, d::Int) = a

function insertentry!(a::Float64, b::AbstractArray, d::Int)
  b[d] = a
  b
end

function evaluategradients!(this::TensorProductMonomialBasis{D,T},
  points::AbstractVector{Point{D}}, v::AbstractArray{TG,2}) where {D,T,TG}
  tpcoor = i -> [ Point{1}(p[i]) for p in points]#@fverdugo [...] allocates a temporary array
  cooruv = [tpcoor(i) for i in 1:D]#@fverdugo [...] allocates a temporary array
  univals = [evaluate(this.univariatebases[i],cooruv[i]) for i in 1:D]#@fverdugo [...] allocates a temporary array
  dervals = [evaluategradients(this.univariatebases[i],cooruv[i]) for i in 1:D]#@fverdugo [...] allocates a temporary array
  cid = ntuple(i -> 1:length(this.univariatebases[i]), D)
  lent = length(T)
  cid = (cid..., 1:lent)
  cid = CartesianIndices(cid)
  E = eltype(T)
  aux = zero(MVector{D,E})
  eb = zero(mutable(T))
  for (i,I) in enumerate(cid)
    d = I[D+1]
    for (p,P) in enumerate(points)
      tpder!(aux, I, p, univals, dervals)
      # @santiagobadia : Any better solution?
      eb = insertentry!(one(E),eb,d)
      v[i,p] = outer(aux,eb)
      insertentry!(zero(E),aux,d)
      insertentry!(zero(E),eb,d)
    end
  end
end

function tpder!(aux::MVector{D,E}, I::CartesianIndex{L}, p::Int, univals, dervals) where {D,E,L}
  for α in 1:D
    val = one(E)
    for β in 1:L-1
      if β != α
        val = val*univals[β][I[β],p][1]
      else
        val = val*dervals[β][I[β],p][1]
      end
    end
    aux[α] = val
  end
end

struct PolynomialBasisWithChangeOfBasis{D,T} <: PolynomialBasis{D,T}
  basis::PolynomialBasis{D,T}
  changeofbasis::Array{Float64,2}
end

function Base.length(this::PolynomialBasisWithChangeOfBasis{D,T})::Int where {D,T}
  length(this.basis)
end

function evaluate!(this::PolynomialBasisWithChangeOfBasis{D,T},
  points::AbstractVector{Point{D}}, v::AbstractArray{T,2}) where {D,T}
  evaluate!(this.basis,points,v)
  v .= this.changeofbasis*v
end

function evaluategradients!(this::PolynomialBasisWithChangeOfBasis{D,T},
  points::AbstractVector{Point{D}}, v::AbstractArray{TG,2}) where {D,T,TG}
  evaluategradients!(this.basis,points,v)
  v .= this.changeofbasis*v
end

end # module Polynomials
