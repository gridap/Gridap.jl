#################################################
# Polynomial family/types with their parameters #
#################################################

"""
Abstract type for 1- or n-D polynomial family of maximum order k.

Implements [`isHierarchical`](@ref).
"""
abstract type Polynomial  <: Field end

"""
Return true if the basis of order `k` of the given `<:Polynomial` type is the union
of the basis of order `k-1` and an other order `k` polynomial.
"""
isHierarchical(::Polynomial) = @abstractmethod

"""
    PolynomialBasis{D,V,K,PT<:Polynomial} <: AbstractVector{PT}

Abstract type representing a generic multivariate polynomial basis.
The parameters are:
- `D`: the spatial dimension
- `V`: the image values type, of type `<:Real` or `<:MultiValue`
- `K`: the maximum order of a basis polynomial in a spatial component
- `PT <: Polynomial`: the polynomial family (must be a concrete type)
"""
abstract type PolynomialBasis{D,V,K,PT<:Polynomial} <: AbstractVector{PT}  end

@inline Base.size(a::PolynomialBasis{D,V}) where {D,V} = @abstractmethod
@inline Base.getindex(a::PolynomialBasis{D,V,K,PT}, i::Integer) where {D,V,K,PT} = PT()
@inline Base.IndexStyle(::PolynomialBasis) = IndexLinear()
@inline return_type(::PolynomialBasis{D,V}) where {D,V} = V

"""
    get_order(b::PolynomialBasis{D,V,K) = K

Return the maximum polynomial order in a dimension, or `0` in 0D.
"""
get_order(::PolynomialBasis{D,V,K}) where {D,V,K} = K


###############################
# 1D internal polynomial APIs #
###############################

"""
    _evaluate_1d!(PT::Type{<:Polynomial},::Val{K},v,x,d)

Evaluates in place the 1D basis polynomials of the given type at one nD point `x`
along the given coordinate 1 ≤ `d` ≤ nD.

`v` is an AbstractMatrix of size (at least} `d`×(`K`+1), such that the 1 ≤ i ≤ `k`+1
values are stored in `v[d,i]`.
"""
function _evaluate_1d!(::Type{<:Polynomial},::Val{K},v::AbstractMatrix{T},x,d) where {K,T<:Number}
  @abstractmethod
end

"""
    _gradient_1d!(PT::Type{<:Polynomial},::Val{K},g,x,d)

Like [`_evaluate_1d!`](@ref), but computes the first derivative of the basis functions.
"""
function _gradient_1d!(::Type{<:Polynomial},::Val{K},g::AbstractMatrix{T},x,d) where {K,T<:Number}
  @abstractmethod
end

"""
    _hessian_1d!(PT::Type{<:Polynomial},::Val{K},g,x,d)

Like [`_evaluate_1d!`](@ref), but computes the second derivative of the basis functions.
"""
function _hessian_1d!(::Type{<:Polynomial},::Val{K},h::AbstractMatrix{T},x,d) where {K,T<:Number}
  @abstractmethod
end

"""
    _derivatives_1d!(PT::Type{<:Polynomial}, ::Val{K}, (v,g,...), x, d)

Same as calling
```
_evaluate_1d!(PT, Val(K), v, x d)
_gradient_1d!(PT, Val(K), g, x d)
          ⋮
```
but with possible performance optimization.
"""
function _derivatives_1d!(  ::Type{<:Polynomial},::Val{K},t::NTuple{N},x,d) where {K,N}
  @abstractmethod
end

function _derivatives_1d!(PT::Type{<:Polynomial},::Val{K},t::NTuple{1},x,d) where K
  _evaluate_1d!(PT, Val(K), t[1], x, d)
end

function _derivatives_1d!(PT::Type{<:Polynomial},::Val{K},t::NTuple{2},x,d) where K
  _evaluate_1d!(PT, Val(K), t[1], x, d)
  _gradient_1d!(PT, Val(K), t[2], x, d)
end

function _derivatives_1d!(PT::Type{<:Polynomial},::Val{K},t::NTuple{3},x,d) where K
  _evaluate_1d!(PT, Val(K), t[1], x, d)
  _gradient_1d!(PT, Val(K), t[2], x, d)
  _hessian_1d!( PT, Val(K), t[3], x, d)
end

# Optimizing evaluation at a single point

function return_cache(f::AbstractVector{PT},x::Point) where PT<:Polynomial
  xs = [x]
  cf = return_cache(f,xs)
  v = evaluate!(cf,f,xs)
  r = CachedArray(zeros(eltype(v),(size(v,2),)))
  r, cf, xs
end

function evaluate!(cache,f::AbstractVector{PT},x::Point) where PT<:Polynomial
  r, cf, xs = cache
  xs[1] = x
  v = evaluate!(cf,f,xs)
  ndof = size(v,2)
  setsize!(r,(ndof,))
  a = r.array
  copyto!(a,v)
  a
end

function return_cache(
  f::FieldGradientArray{N,<:AbstractVector{PT}}, x::Point) where {N,PT<:Polynomial}

  xs = [x]
  cf = return_cache(f,xs)
  v = evaluate!(cf,f,xs)
  r = CachedArray(zeros(eltype(v),(size(v,2),)))
  r, cf, xs
end

function evaluate!(
  cache, f::FieldGradientArray{N,<:AbstractVector{PT}}, x::Point) where {N,PT<:Polynomial}

  r, cf, xs = cache
  xs[1] = x
  v = evaluate!(cf,f,xs)
  ndof = size(v,2)
  setsize!(r,(ndof,))
  a = r.array
  copyto!(a,v)
  a
end

