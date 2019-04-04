module Polynomials

using Numa #@fverdugo to be eliminated
using Numa.Helpers
using Numa.FieldValues
# using Numa.Quadratures
# using Base.Cartesian

export MultivariatePolynomialBasis
export TensorProductMonomialBasis
export MPB_WithChangeOfBasis
# export UnivariatePolynomialBasis
# export UnivariateMonomialBasis

export evaluate
export gradient, ∇
export evaluate! # @santiagobadia : I would not export it
# @fverdugo, basically the clients of polynomial will call evaluate!, gradient, and length

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

# @fverdugo the function we want to overwrite is evaluate! not evaluate
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
 error("Gradient of a gradient not available")
end

const ∇ = gradient

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
Compute the number-th derivative of a monomial at a set of 1D points,
returning an array with first axis basis function label, second axis point label
"""
function derivative(this::UnivariateMonomialBasis,
  points::AbstractVector{Point{1}}; numd=1::Int)::Array{VectorValue{1},2}
  dbas = length(this)
  #@fverdugo we are allocating a vector, thus it is not efficient to call this function in a loop
  # (in unfitted methods it is typically needed to evaluate polynomials within a loop over cut cells)
  v = Array{VectorValue{1},2}(undef, dbas, length(points))
  for (j,p) ∈ enumerate(points)
    for i in 1:length(this)
      val = (i<=numd) ? 0.0 : prod([i-k-1 for k=0:numd-1])p[1]^(i-numd-1)#@fverdugo [...] allocates a temporary array
      v[i,j] = VectorValue{1}(val)
    end
  end
  return v
end

# @fverdugo code repetition in derivative and evaluategradients that can be easily fixed
function evaluategradients!(this::UnivariateMonomialBasis,
  points::AbstractVector{Point{1}},v::AbstractArray{VectorValue{1},2})
  numd = 1 # Changing this number, we get higher order derivatives, to be used
	# for Hessians, etc.
  for (j,p) ∈ enumerate(points)
    for i in 1:length(this)
      val = (i<=numd) ? 0.0 : prod([i-k-1 for k=0:numd-1])p[1]^(i-numd-1)#@fverdugo [...] allocates a temporary array
      v[i,j] = VectorValue{1}(val)
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


Base.length(::Type{ScalarValue}) = 1 #@fverdugo why it is needed? length(::Type{Float64}) already defined by julia
# @santiagobadia : If I comment this line, it does not work
function Base.length(this::TensorProductMonomialBasis{D,T})::Int where {D,T}
  length(T)*prod([length(this.univariatebases[i]) for i in 1:D])#@fverdugo [...] allocates a temporary array
end

function evaluate!(this::TensorProductMonomialBasis{D,T},
  points::AbstractVector{Point{D}}, v::AbstractArray{T,2}) where {D,T}
  tpcoor = i -> [ Point{1}(p[i]) for p in points]#@fverdugo [...] allocates a temporary array
  cooruv = [tpcoor(i) for i in 1:D]#@fverdugo [...] allocates a temporary array
  univals = [evaluate(this.univariatebases[i],cooruv[i]) for i in 1:D]#@fverdugo [...] allocates a temporary array
	cid = ntuple(i -> 1:length(this.univariatebases[i]), D)
	lent = length(T)
	cid = (cid..., 1:lent)
	cid = CartesianIndices(cid)
	for (i,j) in enumerate(cid)
		d = j[D+1]
		for k in 1:length(points)
			val = prod([ univals[i][j[i],k] for i in 1:D ])#@fverdugo [...] allocates a temporary array
			v[i,k] = T(ntuple(i->(i==d) ? val : 0.0, lent)...)# @fverdugo I would say that this is not type stable
		end
	end
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
  for (i,I) in enumerate(cid)
    d = I[D+1]
    for (p,P) in enumerate(points)
			aux = VectorValue{D}([tpder(α, I, p, univals, dervals) for α in 1:D])
			eb = T(ntuple(i->(i==d) ? 1 : 0.0, lent)...)# @fverdugo I would say that this is not type stable
			v[i,p] = outer(aux,eb)
    end
  end
end

function tpder(α::Int, I::CartesianIndex{D}, p::Int, univals, dervals)::Float64 where D
	aux = 1
	for β in 1:D-1
		if β != α
			aux = aux*univals[β][I[β],p][1]
		else
			aux = aux*dervals[β][I[β],p][1]
		end
	end
	return aux
end

struct MPB_WithChangeOfBasis{D,T} <: MultivariatePolynomialBasis{D,T}
  basis::MultivariatePolynomialBasis{D,T}
	changeofbasis::Array{Float64,2}
end

function Base.length(this::MPB_WithChangeOfBasis{D,T})::Int where {D,T}
  length(this.basis)
end

function evaluate!(this::MPB_WithChangeOfBasis{D,T},
  points::AbstractVector{Point{D}}, v::AbstractArray{T,2}) where {D,T}
	evaluate!(this.basis,points,v)
	println(v)
	v .= this.changeofbasis*v
end

# function evaluate(this::MPB_WithChangeOfBasis{D,T},points::AbstractVector{Point{D}}) where {D,T}
#   vals = Array{T,2}(undef,(length(this),length(points)))
#   evaluate!(this.basis,points,vals)
#   return this.changeofbasis*vals
# end
# @fverdugo delete PolynomialsMethods.jl if not needed

# @santiagobadia : Missing evaluategradients for TensorProductMonomialBasis !!!

end # module Polynomials
