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

#"""
#    get_order(::Type{PolynomialType{K}}) where K = K
#"""
#get_order(::Type{PolynomialType{K}}) where K = K

#"""
#    Jacobi{α, β} <: Polynomial
#
#Type representing Jacobi polynomials of parameters `α` and `β`
#"""
#struct Jacobi{α, β} <: Polynomial end
#isHierarchical(::Jacobi) = true


###############################
# 1D internal polynomial APIs #
###############################

# TODO pass order argument as Val{k} to use compile time dispatch to optimize
# edge cases (k=0,1) and possible polynomial coefficients that can be
# pre-computed at compile time

"""
    _evaluate_1d!(::Type{<:Polynomial}, k v, x, d)

Evaluates in place the 1D basis polynomials of the given type at one nD point `x`
along the given coordinate 1 ≤ `d` ≤ nD.

`v` is an AbstractMatrix of size (at least} d×(k+1), such that the 1 ≤ i ≤ `k`+1
values are stored in `v[d,i]`.
"""
function _evaluate_1d!( ::Type{<:Polynomial},k,v::AbstractMatrix{T},x,d) where T<:Number
  @abstractmethod
end

"""
    _gradient_1d!(::Type{<:Polynomial}, k, g, x, d)

Like [`_evaluate_1d!`](@ref), but computes the first derivative of the basis functions.
"""
function _gradient_1d!( ::Type{<:Polynomial},k,g::AbstractMatrix{T},x,d) where T<:Number
  @abstractmethod
end

"""
    _hessian_1d!(::Type{<:Polynomial}, k, h, x, d)

Like [`_evaluate_1d!`](@ref), but computes the second derivative of the basis functions.
"""
function _hessian_1d!( ::Type{<:Polynomial},k,h::AbstractMatrix{T},x,d) where T<:Number
  @abstractmethod
end

"""
    _derivatives_1d!(PT::Type{<:Polynomial}, k, (v,g,...), x, d)

Same as calling
```
_evaluate_1d!(PT, k, v, x d)
_gradient_1d!(PT, k, g, x d)
          ⋮
```
but with possible performence optimization.
"""
function _derivatives_1d!(PT::Type{<:Polynomial},k,t::NTuple{N},x,d) where N
  @abstractmethod
end

function _derivatives_1d!(PT::Type{<:Polynomial},k,t::NTuple{1},x,d)
  _evaluate_1d!(PT, k, t[1], x, d)
end

function _derivatives_1d!(PT::Type{<:Polynomial},k,t::NTuple{2},x,d)
  _evaluate_1d!(PT, k, t[1], x, d)
  _gradient_1d!(PT, k, t[2], x, d)
end

function _derivatives_1d!(PT::Type{<:Polynomial},k,t::NTuple{3},x,d)
  _evaluate_1d!(PT, k, t[1], x, d)
  _gradient_1d!(PT, k, t[2], x, d)
  _hessian_1d!( PT, k, t[3], x, d)
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

