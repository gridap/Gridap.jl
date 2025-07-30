############################
# Polynomial family/types  #
############################

"""
    Polynomial  <: Field

Abstract type for polynomial bases families/types. It has trait
[`isHierarchical`](@ref).
"""
abstract type Polynomial  <: Field end

"""
    isHierarchical(::Type{Polynomial})::Bool

Return `true` if the 1D basis of order `K` of the given [`Polynomial`](@ref)
basis family is the union of the basis of order `K-1` and an other order `K`
polynomial. Equivalently, if the iᵗʰ basis polynomial is of order i-1.

The currently implemented families are [Monomial](@ref), [Legendre](@ref),
[Chebyshev](@ref), [ModalC0](@ref) and [Bernstein](@ref). Only Bernstein is not
hierarchical.
"""
isHierarchical(::Type{<:Polynomial}) = @abstractmethod

testvalue(::Type{PT}) where PT<:Polynomial = isconcretetype(PT) ? PT() : @abstractmethod

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
# ndof_1d: maximum number of 1D monomial in any spatial dimension

"""
    PolynomialBasis{D,V,PT<:Polynomial} <: AbstractVector{PT}

Abstract type representing a generic multivariate polynomial basis.
The parameters are:
- `D`: the spatial dimension
- `V`: the image values type, a concrete type `<:Real` or `<:MultiValue`
- `PT <: Polynomial`: the family of the basis polynomials (must be a concrete type).

The implementations also stores `K`: the maximum order of a basis polynomial in a spatial component
"""
abstract type PolynomialBasis{D,V,PT<:Polynomial} <: AbstractVector{PT}  end

@inline Base.size(::PolynomialBasis{D,V}) where {D,V} = @abstractmethod
@inline Base.getindex(::PolynomialBasis{D,V,PT}, i::Integer) where {D,V,PT} = PT()
@inline Base.IndexStyle(::PolynomialBasis) = IndexLinear()
@inline return_type(::PolynomialBasis{D,V}) where {D,V} = V

"""
    get_order(b::PolynomialBasis)

Return the maximum polynomial order in a dimension, or `0` in 0D.
"""
@inline get_order(::PolynomialBasis) = @abstractmethod
get_order(f::Fields.LinearCombinationFieldVector) = get_order(f.fields)
get_order(f::AbstractVector{<:ConstantField}) = 0

testvalue(::Type{<:PolynomialBasis}) = @abstractmethod


###########
# Helpers #
###########

_q_filter( e,order)  = (maximum(e,init=0) <= order) # ℚₙ
_qh_filter(e,order)  = (maximum(e,init=0) == order) # ℚ̃ₙ = ℚₙ\ℚ₍ₙ₋₁₎
_p_filter( e,order)  = (sum(e) <= order)            # ℙₙ
_ph_filter(e,order)  = (sum(e) == order)            # ℙ̃ₙ = ℙₙ\ℙ₍ₙ₋₁₎
_ser_filter(e,order) = (sum( [ i for i in e if i>1 ] ) <= order) # Serendipity

function _define_terms(filter,orders)
  t = orders .+ 1
  g = (0 .* orders) .+ 1
  cis = CartesianIndices(t)
  co = CartesianIndex(g)
  maxorder = maximum(orders, init=0)
  [ ci for ci in cis if filter(Int[Tuple(ci-co)...],maxorder) ]
end


#########################################
# Generic array of field implementation #
#########################################

function _return_cache(
  f::PolynomialBasis{D}, x,::Type{G},::Val{N_deriv}) where {D,G,N_deriv}

  T = eltype(G)
  np = length(x)
  ndof = length(f)
  ndof_1d = get_order(f) + 1
  # Cache for the returned array
  r = CachedArray(zeros(G,(np,ndof)))
  # Mutable cache for one N_deriv's derivative of a T-valued scalar polynomial
  s = MArray{Tuple{Vararg{D,N_deriv}},T}(undef)
  # Cache for the 1D basis function values in each dimension (to be
  # tensor-producted), and of their N_deriv'th 1D derivatives
  t = ntuple( _ -> CachedArray(zeros(T,(D,ndof_1d))), Val(N_deriv+1))
  (r, s, t...)
end

function _return_val_eltype(b::PolynomialBasis{D,V}, x::AbstractVector{<:Point}) where {D,V}
  xi = testitem(x)
  zVc = zero(eltype(V))
  zxic = zero(eltype(xi))
  T = typeof(zVc*zxic)
  change_eltype(V, T) # Necessary for dual number probagation for autodiff
end

function return_cache(f::PolynomialBasis{D,V}, x::AbstractVector{<:Point}) where {D,V}
  @assert D == length(eltype(x)) "Incorrect number of point components"
  Vr = _return_val_eltype(f,x)
  _return_cache(f,x,Vr,Val(0))
end

function return_cache(
  fg::FieldGradientArray{N,<:PolynomialBasis{D,V}},
  x::AbstractVector{<:Point}) where {N,D,V}

  @assert D == length(eltype(x)) "Incorrect number of point components"
  f = fg.fa
  xi = testitem(x)
  G = _return_val_eltype(f,x)
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

"""
    _get_static_parameters(::PolynomialBasis)

Return a (tuple of) static parameter(s) appended to low level `[...]_nd!` evaluation
calls, default is `Val(get_order(b))`.
"""
_get_static_parameters(b::PolynomialBasis) = Val(get_order(b))

function evaluate!(cache,
  f::PolynomialBasis,
  x::AbstractVector{<:Point})

  r, _, c = cache
  np = length(x)
  _setsize!(f,np,r,c)
  params = _get_static_parameters(f)
  for i in 1:np
    @inbounds xi = x[i]
    _evaluate_nd!(f,xi,r,i,c,params)
  end
  r.array
end

function evaluate!(cache,
  fg::FieldGradientArray{1,<:PolynomialBasis{D,V}},
  x::AbstractVector{<:Point}) where {D,V}

  f = fg.fa
  r, s, c, g = cache
  np = length(x)
  _setsize!(f,np,r,c,g)
  params = _get_static_parameters(f)
  for i in 1:np
    @inbounds xi = x[i]
    _gradient_nd!(f,xi,r,i,c,g,s,params)
  end
  r.array
end

function evaluate!(cache,
  fg::FieldGradientArray{2,<:PolynomialBasis{D,V}},
  x::AbstractVector{<:Point}) where {D,V}

  f = fg.fa
  r, s, c, g, h = cache
  np = length(x)
  _setsize!(f,np,r,c,g,h)
  params = _get_static_parameters(f)
  for i in 1:np
    @inbounds xi = x[i]
    _hessian_nd!(f,xi,r,i,c,g,h,s,params)
  end
  r.array
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
# nD internal polynomial APIs #
###############################

"""
    _evaluate_nd!(b,xi,r,i,c,params)

Compute and assign: `r`[`i`] = `b`(`xi`) = (`b`₁(`xi`), ..., `b`ₙ(`xi`))

where n = length(`b`) (cardinal of the basis), that is the function computes
the basis polynomials at a single point `xi` and sets the result in the `i`th
row of `r`.

- `c` is an implementation specific cache for temporary computation of `b`(`xi`).
- `params` is an optional (tuple of) parameter(s) returned by [`_get_static_parameters(b)`](@ref _get_static_parameters)
"""
function _evaluate_nd!(b::PolynomialBasis, xi, r::AbstractMatrix, i, c, params)
  @abstractmethod
end

"""
    _gradient_nd!(b,xi,r,i,c,g,s,params)

Compute and assign: `r`[`i`] = ∇`b`(`xi`) = (∇`b`₁(`xi`), ..., ∇`b`ₙ(`xi`))

where n = length(`b`) (cardinal of the basis), like [`_evaluate_nd!`](@ref) but
for gradients of `b`ₖ(`xi`), and

- `g` is an implementation specific cache for temporary computation of `∇b`(`xi`).
- `s` is a mutable length `D` cache for ∇`b`ₖ(`xi`).
"""
function _gradient_nd!(b::PolynomialBasis, xi, r::AbstractMatrix, i, c, g, s::MVector, params)
  @abstractmethod
end

"""
    _hessian_nd!(b,xi,r,i,c,g,h,s,params)

Compute and assign: `r`[`i`] = H`b`(`xi`) = (H`b`₁(`xi`), ..., H`b`ₙ(`xi`))

where n = length(`b`) (cardinal of the basis), like [`_evaluate_nd!`](@ref) but
for hessian matrices/tensor of `b`ₖ(`xi`), and

- `h` is an implementation specific cache for temporary computation of `∇∇b`(`xi`).
- `s` is a mutable `D`×`D` cache for H`b`ₖ(`xi`).
"""
function _hessian_nd!(b::PolynomialBasis, xi, r::AbstractMatrix, i, c, g, h, s::MMatrix, params)
  @abstractmethod
end


###############################
# 1D internal polynomial APIs #
###############################

"""
    _evaluate_1d!(PT::Type{<:Polynomial},::Val{K},c,x,d)

Evaluates in place the 1D basis polynomials of the family `PT` at one D-dim.
point `x` along the given coordinate 1 ≤ `d` ≤ D.

`c` is an AbstractMatrix of size (at least) `d`×(`K`+1), such that the
1 ≤ i ≤ `k`+1 values are stored in `c[d,i]`.
"""
function _evaluate_1d!(::Type{<:Polynomial},::Val{K},c::AbstractMatrix{T},x,d) where {K,T<:Number}
  @abstractmethod
end

"""
    _gradient_1d!(PT::Type{<:Polynomial},::Val{K},g,x,d)

Like [`_evaluate_1d!`](@ref), but computes the first derivative of the basis
polynomials.
"""
function _gradient_1d!(::Type{<:Polynomial},::Val{K},g::AbstractMatrix{T},x,d) where {K,T<:Number}
  @abstractmethod
end

"""
    _hessian_1d!(PT::Type{<:Polynomial},::Val{K},g,x,d)

Like [`_evaluate_1d!`](@ref), but computes the second derivative of the basis
polynomials.
"""
function _hessian_1d!(::Type{<:Polynomial},::Val{K},h::AbstractMatrix{T},x,d) where {K,T<:Number}
  @abstractmethod
end

# Dispatch helpers for base cases
const Val_01  = Union{Val{0},Val{1}}
const Val_012 = Union{Val{0},Val{1},Val{2}}

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
function _derivatives_1d!(  ::Type{<:Polynomial},v::Val,t::NTuple{N},x,d) where N
  @abstractmethod
end

function _derivatives_1d!(PT::Type{<:Polynomial},v::Val,t::NTuple{1},x,d)
  @inline _evaluate_1d!(PT, v, t[1], x, d)
end

function _derivatives_1d!(PT::Type{<:Polynomial},v::Val,t::NTuple{2},x,d)
  @inline _evaluate_1d!(PT, v, t[1], x, d)
  @inline _gradient_1d!(PT, v, t[2], x, d)
end

function _derivatives_1d!(PT::Type{<:Polynomial},v::Val,t::NTuple{3},x,d)
  @inline _evaluate_1d!(PT, v, t[1], x, d)
  @inline _gradient_1d!(PT, v, t[2], x, d)
  @inline _hessian_1d!( PT, v, t[3], x, d)
end
