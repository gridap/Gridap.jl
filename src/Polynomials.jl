module Polynomials

using Numa #@fverdugo to be eliminated
using Numa.Helpers
using Numa.FieldValues
# using Numa.Quadratures
# using Base.Cartesian

export MultivariatePolynomialBasis
# export TensorProductMonomialBasis
# export UnivariatePolynomialBasis
# export UnivariateMonomialBasis

# export evaluate!
export evaluate
export gradient, ∇

"""
Abstract type representing a multivariate polynomial basis
with value of type T in a coordinate space of D dimensions
"""
abstract type MultivariatePolynomialBasis{D,T} end

Base.length(::MultivariatePolynomialBasis)::Int = @abstractmethod


"""
Same as evaluate! but allocates output
"""
function evaluate(self::MultivariatePolynomialBasis{D,T},points::AbstractArray{Point{D},1}) where {D,T}
  vals = Array{T,2}(undef,(length(self),length(points)))
  evaluate!(self,points,vals)
  vals
end

"""
First axis of v for dofs, second for points
"""
evaluate!(::MultivariatePolynomialBasis{D,T},::AbstractVector{Point{D}},v::AbstractArray{T,2}) where {D,T} = @abstractmethod

"""
Object that represents the gradient of the elements of a basis
"""
struct GradMultivariatePolynomialBasis{D,TG,T,B<:MultivariatePolynomialBasis{D,T}} <: MultivariatePolynomialBasis{D,TG}
  basis::B
end

"""
evaluate! overwritten for gradients of bases
"""
function evaluate(self::GradMultivariatePolynomialBasis{D,T},points::AbstractArray{Point{D},1}) where {D,T}
  vals = Array{T,2}(undef,(length(self),length(points)))
  evaluategradients!(self.basis,points,vals)
  vals
end

"""
Create a `GradMultivariatePolynomialBasis` object from a basis
`MultivariatePolynomialBasis`. The result is a
MultivariatePolynomialBasis{TG,D} where TG is a type whose rank is one unit
greater than the one of T
"""
function gradient(self::MultivariatePolynomialBasis{D,T}) where{D,T}
  TG = outer(Point{D},T)
  B = typeof(self)
  GradMultivariatePolynomialBasis{D,TG,T,B}(self)
end

Base.length(this::GradMultivariatePolynomialBasis)::Int = length(this.basis)

function gradient(::GradMultivariatePolynomialBasis{D,T,B} where{D,T,B})
 @error("Gradient of a gradient not available") end

const ∇ = gradient

Base.:*(::typeof(∇),f) = ∇(f)

"""
Abstract type representing a univariate polynomial basis in dimension one
"""
abstract type UnivariatePolynomialBasis <: MultivariatePolynomialBasis{1,ScalarValue} end

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

"""
Compute the numder-th derivative of a monomial at a set of 1D points,
returning an array with first axis basis function label, second axis point label
"""
function derivative(this::UnivariateMonomialBasis,
  points::AbstractVector{Point{1}}; numd=1::Int)::Array{VectorValue{1},2}
  dbas = length(this)
  v = Array{VectorValue{1},2}(undef, dbas, length(points))
  for (j,p) ∈ enumerate(points)
    for i in 1:length(this)
      val = (i<=numd) ? 0.0 : prod([i-k-1 for k=0:numd-1])p[1]^(i-numd-1)
      v[i,j] = VectorValue{1}(val)
    end
  end
  return v
end

function evaluategradients!(this::UnivariateMonomialBasis,
  points::AbstractVector{Point{1}},v::AbstractArray{VectorValue{1},2})
  derivative(this, points)
  numd = 1
  for (j,p) ∈ enumerate(points)
    for i in 1:length(this)
      val = (i<=numd) ? 0.0 : prod([i-k-1 for k=0:numd-1])p[1]^(i-numd-1)
      v[i,j] = VectorValue{1}(val)
    end
  end
end

"""
Multivariate monomial basis obtained as tensor product of univariate polynomial
basis per dimension
"""
struct TensorProductMonomialBasis{D,T} <: MultivariatePolynomialBasis{D,T}
  univariatebases::Vector{UnivariateMonomialBasis}
end

"""
Provide a `TensorProductMonomialBasis` for a vector `order` providing the order per
dimension
"""
function TensorProductMonomialBasis{D,T}(order::Vector{Int64}) where {D,T}
  @assert(length(order) == D)
  uvmbs = [UnivariateMonomialBasis(order[i]) for i=1:length(order)]
  TensorProductMonomialBasis{D,T}(uvmbs)
end

Base.length(::Type{ScalarValue}) = 1
function Base.length(this::TensorProductMonomialBasis{D,T})::Int where {D,T}
  length(T)*prod([length(this.univariatebases[i]) for i in 1:D])
end

function evaluate!(this::TensorProductMonomialBasis{D,T},
  points::AbstractVector{Point{D}}, v::AbstractArray{T,2}) where {D,T}
  tpcoor = i -> [ Point{1}(p[i]) for p in points]
  cooruv = [tpcoor(i) for i in 1:D]
  univals = [evaluate(this.univariatebases[i],cooruv[i]) for i in 1:D]
  cid = ntuple(i -> 1:length(this.univariatebases[i]), D)
  lent = length(T)
  cid = (cid..., 1:lent)
  cid = CartesianIndices(cid)
  for (i,j) in enumerate(cid)
    d = j[D+1]
    for k in 1:length(points)
      val = prod([ univals[i][j[i],k] for i in 1:D ])
      v[i,k] = T(ntuple(i->(i==d) ? val : 0.0, lent))
    end
  end
end

function oldevaluate!(this::TensorProductMonomialBasis{D,T},
  points::AbstractVector{Point{D}}, v::AbstractArray{T,2}) where {D,T}
  tpcoor = i -> [ Point{1}(p[i]) for p in points]
  cooruv = [tpcoor(i) for i in 1:D]
  univals = [evaluate(this.univariatebases[i],cooruv[i]) for i in 1:D]
  cid = ntuple(i -> 1:length(this.univariatebases[i]), D)
  cid = CartesianIndices(cid)
  for (i,j) in enumerate(cid)
    for k in 1:length(points)
      v[i,k] = prod([ univals[i][j[i],k] for i in 1:D ])
    end
  end
end

end # module Polynomials
