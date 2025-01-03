############################
# Polynomial family/types  #
############################

"""
    Polynomial  <: Field

Abstract type for polynomial bases families/types. It has trait
[`isHierarchical`](@ref).

The currently implemented families are [Monomial](@ref), [Legendre](@ref),
[Chebyshev](@ref), [ModalC0](@ref) and [Bernstein](@ref). Only Bernstein is not
hierarchical.
"""
abstract type Polynomial  <: Field end

"""
    isHierarchical(::Type{Polynomial})

Return true if the 1D basis of order `K` of the given [`Polynomial`](@ref) type
is the union of the basis of order `K-1` and an other order `K` polynomial.
Equivalently, if the iᵗʰ basis polynomial is of order i-1.
"""
isHierarchical(::Type{Polynomial}) = @abstractmethod


###########################################
# Polynomial basis abstract type and APIs #
###########################################

# Notations:
#
# D: spatial / input space dimension
# T: scalar type (Float64, ...)
# V: concrete type of image values (T, VectorValue{D,T} etc.)
# G: concrete MultiValue type holding the gradient or hessian of a function of
#      value V, i.e.  gradient_type(V,Point{D})
#                 or  gradient_type(gradient_type(V,p::Point{D}), p)
#
# PT:      a concrete `Polynomial` type
# K:       integer polynomial order (maximum order of any component and in any direction in nD).
# np:      number of points at which a basis is evaluated
# ndof:    number of basis polynomials
# ndof_1d: maximum of 1D polynomial vector in any spatial dimension

"""
    PolynomialBasis{D,V,K,PT<:Polynomial} <: AbstractVector{PT}

Abstract type representing a generic multivariate polynomial basis.
The parameters are:
- `D`: the spatial dimension
- `V`: the image values type, of type `<:Real` or `<:MultiValue`
- `K`: the maximum order of a basis polynomial in a spatial component
- `PT <: Polynomial`: the polynomial family (must be a concrete type).
"""
abstract type PolynomialBasis{D,V,K,PT<:Polynomial} <: AbstractVector{PT}  end

@inline Base.size(::PolynomialBasis{D,V}) where {D,V} = @abstractmethod
@inline Base.getindex(::PolynomialBasis{D,V,K,PT}, i::Integer) where {D,V,K,PT} = PT()
@inline Base.IndexStyle(::PolynomialBasis) = IndexLinear()
@inline return_type(::PolynomialBasis{D,V}) where {D,V} = V

"""
    get_order(b::PolynomialBasis{D,V,K}) = K

Return the maximum polynomial order in a dimension, is should be `0` in 0D.
"""
@inline get_order(::PolynomialBasis{D,V,K}) where {D,V,K} = K


################################
# Generic field implementation #
################################

function _return_cache(
  f::PolynomialBasis{D}, x,::Type{G},::Val{N_deriv}) where {D,G,N_deriv}

  @assert D == length(eltype(x)) "Incorrect number of point components"
  T = eltype(G)
  np = length(x)
  ndof = length(f)
  ndof_1d = get_order(f) + 1
  # Cache for the returned array
  r = CachedArray(zeros(G,(np,ndof)))
  # Cache for the 1D basis function values in each dimension (to be
  # tensor-producted), and of their N_deriv'th 1D derivatives
  t = ntuple( _ -> CachedArray(zeros(T,(D,ndof_1d ))), Val(N_deriv+1))
  (r, t...)
end

function return_cache(f::PolynomialBasis{D,V}, x::AbstractVector{<:Point}) where {D,V}
  _return_cache(f,x,V,Val(0))
end

function return_cache(
  fg::FieldGradientArray{N,<:PolynomialBasis{D,V}},
  x::AbstractVector{<:Point}) where {N,D,V}

  f = fg.fa
  xi = testitem(x)
  G = V
  for _ in 1:N
    G = gradient_type(G,xi)
  end
  _return_cache(f,x,G,Val(N))
end


function _setsize!(f::PolynomialBasis{D}, np, r, t...) where D
  ndof = length(f)
  ndof_1d = get_order(f) + 1
  setsize!(r,(np,ndof))
  for c in t
    setsize!(c,(D,ndof_1d))
  end
end

function evaluate!(cache,
  f::PolynomialBasis,
  x::AbstractVector{<:Point})

  r, c = cache
  np = length(x)
  _setsize!(f,np,r,c)
  for i in 1:np
    @inbounds xi = x[i]
    _evaluate_nd!(f,xi,r,i,c)
  end
  r.array
end

function evaluate!(cache,
  fg::FieldGradientArray{1,<:PolynomialBasis{D,V}},
  x::AbstractVector{<:Point}) where {D,V}

  f = fg.fa
  r, c, g = cache
  np = length(x)
  _setsize!(f,np,r,c,g)
  s = zero(Mutable(VectorValue{D,eltype(V)}))
  for i in 1:np
    @inbounds xi = x[i]
    _gradient_nd!(f,xi,r,i,c,g,s)
  end
  r.array
end

function evaluate!(cache,
  fg::FieldGradientArray{2,<:PolynomialBasis{D,V}},
  x::AbstractVector{<:Point}) where {D,V}

  f = fg.fa
  r, c, g, h = cache
  np = length(x)
  _setsize!(f,np,r,c,g,h)
  s = zero(Mutable(TensorValue{D,D,eltype(V)}))
  for i in 1:np
    @inbounds xi = x[i]
    _hessian_nd!(f,xi,r,i,c,g,h,s)
  end
  r.array
end


###############################
# nD internal polynomial APIs #
###############################

"""
    _evaluate_nd!(b,xi,r,i,c)

Compute and assign: `r`[`i`] = `b`(`xi`) = (`b`₁(`xi`), ..., `b`ₙ(`xi`))

where n = length(`b`) (cardinal of the basis), that is the function computes
the basis polynomials at a single point `xi` and setting the result in the `i`th
row of `r`.
"""
function _evaluate_nd!(
  b::PolynomialBasis, xi,
  r::AbstractMatrix, i,
  c::AbstractMatrix)

  @abstractmethod
end

"""
    _gradient_nd!(b,xi,r,i,c,g,s)

Compute and assign: `r`[`i`] = ∇`b`(`xi`) = (∇`b`₁(`xi`), ..., ∇`b`ₙ(`xi`))

where n = length(`b`) (cardinal of the basis), like [`_evaluate_nd!`](@ref) but
for gradients of `b`ₖ(`xi`), and

- `g` is a mutable `D`×`K` cache (for 1D poly deriv evals).
- `s` is a mutable length `D` cache for ∇`b`ₖ(`xi`).
"""
function _gradient_nd!(
  b::PolynomialBasis, xi,
  r::AbstractMatrix, i,
  c::AbstractMatrix,
  g::AbstractMatrix,
  s::MVector)

  @abstractmethod
end


"""
    _hessian_nd!(b,xi,r,i,c,g,h,s)

Compute and assign: `r`[`i`] = H`b`(`xi`) = (H`b`₁(`xi`), ..., H`b`ₙ(`xi`))

where n = length(`b`) (cardinal of the basis), like [`_evaluate_nd!`](@ref) but
for hessian matrices/tensor of `b`ₖ(`xi`), and

- `h` is a mutable `D`×`K` cache (for 1D poly second deriv evals).
- `s` is a mutable `D`×`D` cache for H`b`ₖ(`xi`).
"""
function _hessian_nd!(
  b::PolynomialBasis, xi,
  r::AbstractMatrix, i,
  c::AbstractMatrix,
  g::AbstractMatrix,
  h::AbstractMatrix,
  s::MMatrix)

  @abstractmethod
end


##############################################
# Optimizing of evaluation at a single point #
##############################################

function return_cache(f::PolynomialBasis,x::Point)
  xs = [x]
  cf = return_cache(f,xs)
  v = evaluate!(cf,f,xs)
  r = CachedArray(zeros(eltype(v),(size(v,2),)))
  r, cf, xs
end

function evaluate!(cache,f::PolynomialBasis,x::Point)
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
  f::FieldGradientArray{N,<:PolynomialBasis}, x::Point) where N

  xs = [x]
  cf = return_cache(f,xs)
  v = evaluate!(cf,f,xs)
  r = CachedArray(zeros(eltype(v),(size(v,2),)))
  r, cf, xs
end

function evaluate!(
  cache, f::FieldGradientArray{N,<:PolynomialBasis}, x::Point) where N

  r, cf, xs = cache
  xs[1] = x
  v = evaluate!(cf,f,xs)
  ndof = size(v,2)
  setsize!(r,(ndof,))
  a = r.array
  copyto!(a,v)
  a
end


###############################
# 1D internal polynomial APIs #
###############################

"""
    _evaluate_1d!(PT::Type{<:Polynomial},::Val{K},c,x,d)

Evaluates in place the 1D basis polynomials of the given type at one nD point `x`
along the given coordinate 1 ≤ `d` ≤ nD.

`c` is an AbstractMatrix of size (at least) `d`×(`K`+1), such that the 1 ≤ i ≤ `k`+1
values are stored in `c[d,i]`.
"""
function _evaluate_1d!(::Type{<:Polynomial},::Val{K},c::AbstractMatrix{T},x,d) where {K,T<:Number}
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
    _derivatives_1d!(PT::Type{<:Polynomial}, ::Val{K}, (c,g,...), x, d)

Same as calling
```
_evaluate_1d!(PT, Val(K), c, x d)
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

