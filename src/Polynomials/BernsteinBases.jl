"""
    Bernstein <: Polynomial

Type representing Bernstein polynomials, c.f. [Bernstein polynomials](@ref) section.
"""
struct Bernstein <: Polynomial end

isHierarchical(::Type{Bernstein}) = false


#####################################
# Cartesian product Bernstein bases #
#####################################

"""
    BernsteinBasis{D,V} = CartProdPolyBasis{D,V,Bernstein}

Alias for cartesian product of a scalar tensor 1D-Bernstein basis, scalar valued
or multivalued.
"""
const BernsteinBasis{D,V} = CartProdPolyBasis{D,V,Bernstein}

"""
    BernsteinBasis(::Val{D}, ::Type{V}, order::Int, terms::Vector)
    BernsteinBasis(::Val{D}, ::Type{V}, order::Int [, filter::Function])
    BernsteinBasis(::Val{D}, ::Type{V}, orders::Tuple [, filter::Function])

High level constructors of [`BernsteinBasis`](@ref).
"""
BernsteinBasis(args...) = CartProdPolyBasis(Bernstein, args...)


################################
# 1D evaluation implementation #
################################

@inline function _de_Casteljau_step_1D!(v,d,i,λ1,λ2)
  # i = k+1

  # vₖ <- xvₖ₋₁            # Bᵏₖ(x) = x*Bᵏ⁻¹ₖ₋₁(x)
  v[d,i] = λ2*v[d,i-1]
  # vⱼ <- xvⱼ₋₁ + (1-x)vⱼ  # Bᵏⱼ(x) = x*Bᵏ⁻¹ⱼ₋₁(x) + (1-x)*Bᵏ⁻¹ⱼ(x) for j = k-1, k-2, ..., 1
  for l in i-1:-1:2
    v[d,l] = λ2*v[d,l-1] + λ1*v[d,l]
  end
  # v₀ <- (1-x)v₀          # Bᵏ₀(x) = (1-x)*Bᵏ⁻¹₀(x)
  v[d,1] = λ1*v[d,1]
  nothing
end

# jth Bernstein poly of order K at x:
# Bᵏⱼ(x) = binom(K,j) * x^j * (1-x)^(K-j) = x*Bᵏ⁻¹ⱼ₋₁(x) + (1-x)*Bᵏ⁻¹ⱼ(x)
function _evaluate_1d!(::Type{Bernstein},K::Int,v::AbstractMatrix{T},x,d) where T<:Number
  @inbounds begin
    if iszero(K)
      v[d,1] = one(T)
      return
    end

    n = K + 1 # n > 1
    λ2 = x[d]
    λ1 = one(T) - λ2

    # In place De Casteljau: init with B¹₀(x)=x and B¹₁(x)=1-x
    v[d,1] = λ1
    v[d,2] = λ2

    for i in 3:n
      _de_Casteljau_step_1D!(v,d,i,λ1,λ2)
    end
  end
  # still optimisable for K > 2/3:
  # - compute bj = binomials(k,j) ∀j at compile time
  # - compute vj = xʲ*(1-x)ᴷ⁻ʲ recursively in place like De Casteljau (saving half the redundant multiplications)
  # - do it in a stack allocated cache (MVector, Bumber.jl)
  # - @simd affect bj * vj in v[d,i] for all j
end

# First derivative of the jth Bernstein poly of order K at x:
# (Bᵏⱼ)'(x) = K * ( Bᵏ⁻¹ⱼ₋₁(x) - Bᵏ⁻¹ⱼ(x) )
#           = K * x^(j-1) * (1-x)^(K-j-1) * ((1-x)*binom(K-1,j-1) - x*binom(K-1,j))
function _gradient_1d!(::Type{Bernstein},K::Int,g::AbstractMatrix{T},x,d) where T<:Number
  @inbounds begin

    if K<2 # base cases
      z = zero(T)
      o = one(T)
      K==0 && (g[d,1] = z)
      K==1 && (g[d,1] =-o;  g[d,2] = o)
      return
    end

    n = K + 1 # n > 2

    # De Casteljau for Bᵏ⁻¹ⱼ for j = k-1, k-2, ..., 1
    _evaluate_1d!(Bernstein,K-1,g,x,d)

    # gₖ <- K*gₖ₋₁         # ∂ₓBᵏₖ(x) = K*Bᵏ⁻¹ₖ₋₁(x)
    g[d,n] = K*g[d,n-1]
    # gⱼ <- K(gⱼ₋₁ + gⱼ)   # ∂ₓBᵏⱼ(x) = K(Bᵏ⁻¹ⱼ₋₁(x) - Bᵏ⁻¹ⱼ(x)) for j = k-1, k-2, ..., 1
    for l in n-1:-1:2
      g[d,l] = K*(g[d,l-1] - g[d,l])
    end
    # g₀ <- K*g₀           # ∂ₓBᵏ₀(x) = -K*Bᵏ⁻¹₀(x)
    g[d,1] = -K*g[d,1]
  end
end


# Second derivative of the jth Bernstein poly of order K at x:
# (Bᵏⱼ)''(x) = K(K-1) * ( Bᵏ⁻²ⱼ₋₂(x) -2*Bᵏ⁻²ⱼ₋₁(x) + Bᵏ⁻²ⱼ(x) )
#            = K(K-1) * x^(j-2) * (1-x)^(K-j-2) * ( (1-x)^2*binom(K-2,j-2)
#                  - 2x*(1-x)*binom(K-2,j-1) + (x)^2*binom(K-2,j)
#              )
function _hessian_1d!(::Type{Bernstein},K::Int,h::AbstractMatrix{T},x,d) where T<:Number
  @inbounds begin

    if K<3 # base cases
      z = zero(T)
      o = one(T)
      K==0 && (h[d,1] = z)
      K==1 && (h[d,1] = z;  h[d,2] = z)
      K==2 && (h[d,1] = 2o; h[d,2] = -4o; h[d,3] =  2o)
      return
    end

    n = K + 1 # n > 3
    KK = K*(K-1)

    # De Casteljau for Bᵏ⁻²ⱼ for j = k-2, k-3, ..., 1
    _evaluate_1d!(Bernstein,K-2,h,x,d)

    # hₖ   <- K(K-1)*hₖ₋₂
    h[d,n] = KK*h[d,n-2]
    # hₖ₋₁ <- K(K-1)*(-2*hₖ₋₁ + hₖ₋₂)
    h[d,n-1] = KK*( h[d,n-3] -2*h[d,n-2] )

    # hⱼ <- K(K-1)(hⱼ₋₂ -2hⱼ₋₁ + hⱼ)
    for l in n-2:-1:3
      h[d,l] = KK*( h[d,l-2] -2*h[d,l-1] + h[d,l] )
    end

    # h₁ <- K(K-1)*(-2h₀ + h₁)
    h[d,2] = KK*( -2*h[d,1] + h[d,2] )
    # h₀ <- K(K-1)*h₀
    h[d,1] = KK*h[d,1]
  end
end

function _derivatives_1d!(::Type{Bernstein},K,t::NTuple{2},x,d)
  if K < 2
    @inline _evaluate_1d!(Bernstein, K, t[1], x, d)
    @inline _gradient_1d!(Bernstein, K, t[2], x, d)
    return
  end

  @inbounds begin
    n = K + 1 # n > 2
    v, g = t

    λ2 = x[d]
    λ1 = one(eltype(v)) - λ2

    # De Casteljau for Bᵏ⁻¹ⱼ for j = k-1, k-2, ..., 1
    _evaluate_1d!(Bernstein,K-1,v,x,d)

    # Compute gradients as _gradient_1d!
    g[d,n] = K*v[d,n-1]
    @simd for l in n-1:-1:2
      g[d,l] = K*(v[d,l-1] - v[d,l])
    end
    g[d,1] = -K*v[d,1]

    # Last step of De Casteljau for _evaluate_1d!
    _de_Casteljau_step_1D!(v,d,n,λ1,λ2)
  end
end

function _derivatives_1d!(::Type{Bernstein},K,t::NTuple{3},x,d)
  if K < 3
    @inline _evaluate_1d!(Bernstein, K, t[1], x, d)
    @inline _gradient_1d!(Bernstein, K, t[2], x, d)
    @inline _hessian_1d!( Bernstein, K, t[3], x, d)
    return
  end

  @inbounds begin
    n = K + 1 # n > 3
    v, g, h = t

    λ2 = x[d]
    λ1 = one(eltype(v)) - λ2

    # De Casteljau until Bᵏ⁻²ⱼ ∀j
    _evaluate_1d!(Bernstein,K-2,v,x,d)

    # Compute hessians as in _hessian_1d!
    KK = K*(K-1)
    h[d,n] = KK*v[d,n-2]
    h[d,n-1] = KK*( v[d,n-3] -2*v[d,n-2] )
    @simd for l in n-2:-1:3
      h[d,l] = KK*( v[d,l-2] -2*v[d,l-1] + v[d,l] )
    end
    h[d,2] = KK*( -2*v[d,1] + v[d,2] )
    h[d,1] = KK*v[d,1]

    # One step of De Casteljau to get Bᵏ⁻¹ⱼ ∀j
    _de_Casteljau_step_1D!(v,d,n-1,λ1,λ2)

    # Compute gradients as in _gradient_1d!
    g[d,n] = K*v[d,n-1]
    @simd for l in n-1:-1:2
      g[d,l] = K*(v[d,l-1] - v[d,l])
    end
    g[d,1] = -K*v[d,1]

    # Last step of De Casteljau for _evaluate_1d!
    _de_Casteljau_step_1D!(v,d,n,λ1,λ2)
  end
end


###################################################
# Bernstein bases on simplices using de Casteljau #
###################################################

"""
    BernsteinBasisOnSimplex{D,V,M} <: PolynomialBasis{D,V,Bernstein}

Type for the multivariate Bernstein basis in barycentric coordinates,
c.f. [Bernstein polynomials](@ref) section of the documentation. If `V` is not
scalar, a Cartesian product of the Bernstein scalar basis is made for each
independent component of `V`.

The index of `B_α` in the basis is [`bernstein_term_id(α)`](@ref bernstein_term_id).

`M` is Nothing for the reference tetrahedra barycentric coordinates or
`SMatrix{D+1,D+1}` if some simplex (triangle, tetrahedra, ...) vertices
coordinates are given.
"""
struct BernsteinBasisOnSimplex{D,V,M,K} <: PolynomialBasis{D,V,Bernstein}
  cart_to_bary_matrix::M #  Nothing or SMatrix{D+1,D+1}

  function BernsteinBasisOnSimplex{D}(::Type{V},order::Int,vertices=nothing) where {D,V}
    _simplex_vertices_checks(Val(D), vertices)

    cart_to_bary_matrix = _compute_cart_to_bary_matrix(vertices, Val(D+1))
    M = typeof(cart_to_bary_matrix) # Nothing or SMatrix
    K = order
    new{D,V,M,K}(cart_to_bary_matrix)
  end
end

function _simplex_vertices_checks(::Val{D}, vertices) where D
  if !isnothing(vertices)
    @check length(vertices) == D+1 "$D+1 vertices are required to define a $D-dim simplex, got $(length(vertices))"
    @check eltype(vertices) <: Point{D} "Vertices should be of type <:Point{$D}, got $(eltype(vertices))"
  end
end

"""
    BernsteinBasisOnSimplex(::Val{D},::Type{V},order::Int)
    BernsteinBasisOnSimplex(::Val{D},::Type{V},order::Int,vertices)

Constructors for [`BernsteinBasisOnSimplex`](@ref).

If specified, `vertices` is a collection of `D+1` `Point{D}` defining a simplex
used to compute the barycentric coordinates from, it must be non-degenerated
(have nonzero volume).
"""
function BernsteinBasisOnSimplex(::Val{D},::Type{V},order::Int,vertices=nothing) where {D,V}
  BernsteinBasisOnSimplex{D}(V,order,vertices)
end

Base.size(b::BernsteinBasisOnSimplex{D,V}) where {D,V} = (num_indep_components(V)*binomial(D+get_order(b),D),)
get_order(::BernsteinBasisOnSimplex{D,V,M,K}) where {D,V,M,K} = K
get_orders(b::BernsteinBasisOnSimplex{D}) where D = tfill(get_order(b), Val(D))

_get_parameters(::BernsteinBasisOnSimplex{D,V,M,K}) where {D,V,M,K} = Val(K)

function testvalue(::Type{BernsteinBasisOnSimplex{D,V,M,K}}) where {D,V,M,K}
  if M == Nothing
    vertices = nothing
  else
    Pt = Point{D,eltype(M)}
    vertices = ntuple( j -> Pt( ntuple( i -> j==i+1, Val(D)) ), Val(D+1))
  end
  BernsteinBasisOnSimplex{D}(V,K,vertices)
end


#####################
# Bernstein Helpers #
#####################

"""
    _compute_cart_to_bary_matrix(vertices, ::Val{N})
    _compute_cart_to_bary_matrix(::Nothing,::Val) = nothing

For the given the vertices of a `D`-simplex (`D` = `N`-1), computes the change
of coordinate matrix `x_to_λ` from cartesian to barycentric, such that
`λ` = `x_to_λ` * `x` with `sum(λ) == 1` and `x == sum(λ .* vertices)`.
"""
function _compute_cart_to_bary_matrix(vertices, ::Val{N}) where N
  T = eltype(eltype(vertices))
  λ_to_x = MMatrix{N,N,T}(undef)
  for (i,v) in enumerate(vertices)
    λ_to_x[:,i] .= tuple(one(T), v...)
  end

  x_to_λ = inv(λ_to_x) # Doesn't throw if singular because this is a StaticArrays Matrix
  msg =  "The simplex defined by the given vertices is degenerated (is flat / has zero volume)."
  !all(isfinite, x_to_λ) && throw(DomainError(vertices,msg))

  return SMatrix{N,N,T}(x_to_λ)
end
_compute_cart_to_bary_matrix(::Nothing, ::Val) = nothing

"""
    _cart_to_bary(x::Point{D,T}, ::Nothing)

Compute the barycentric coordinates with respect to the reference simplex of the
given cartesian coordinates `x`, that is `λ`=(x1, ..., xD, 1-x1-x2-...-xD).
"""
@inline function _cart_to_bary(x::Point{D,T}, ::Nothing) where {D,T}
  sum_x = sum(x,init=zero(T))
  return SVector(1-sum_x, x...)
end

"""
    _cart_to_bary(x::Point{D,T}, x_to_λ)

Compute the barycentric coordinates of the given cartesian coordinates `x` using
the `x_to_λ` change of coordinate matrix, see [`_compute_cart_to_bary_matrix`](@ref).
"""
@inline function _cart_to_bary(x::Point{D,T}, x_to_λ) where {D,T}
  x_1 = SVector{D+1,T}(one(T), x...)
  return x_to_λ*x_1
end

"""
    bernstein_terms(K,D)

Return the vector of multi-indices for the `D`-dimensional Bernstein basis of
order `K`, that is

``Iₖᴰ = { α ∈ {0:K}ᴰ⁺¹ | |α| = K }``

ordered in decreasing lexicographic order, e.g. {200, 110, 101, 020, 011, 002}
for K=2, D=2.
"""
function bernstein_terms(K,D)
  terms = collect(multiexponents(D+1,K))
  terms = convert(Vector{Vector{Int}}, terms)
end


########################
# de Casteljau helpers #
########################

"""
    bernstein_term_id(α)

For a given Bernstein multi-index `α` (vector or tuple), return the associated
linear index of `α` ordered in decreasing lexicographic order, that is the `i`
such that

    (i,α) ∈ enumerate(bernstein_terms(K,D))

where K = sum(`α`), see also [`bernstein_terms`](@ref).
"""
function bernstein_term_id(α)
  D = length(α)-1
  @check D ≥ 0
  @inbounds i = sum(_L_slices_size(L, D, _L_slice_2(L,α,D)) for L in 1:D; init=0) + 1
  return i
end


# For a given Bernstein term `α`, return the index (starting from 0) of the
# (D-`L`)-slice to which `α` belongs within the (D-`L`-1)-slice of
# the D-multiexponent simplex (D = `N`-1).
#
# In a D-multiexponent simplex of elements `α`, ordered in a vector in
# decreasing lexicographic order, the (D-`L`)-slices are the consecutive `α`s
# having iddentical first `L` indices `α`.
"""
    _L_slice(L, α, D)

where `L` ∈ 1:N, returns `sum( α[i] for i in (L+1):(D+1) )`.
"""
Base.@propagate_inbounds _L_slice_2(L, α, D) = sum( α[i] for i in (L+1):(D+1) )

# Return the length of the `l`-1 first (`D`-`L`)-slices in the
# `D`-multiexponent simplex (ordered in decreasing lexicographic order).
# Those numbers are the "(`D`-`L`)-simplex numbers".
"""
    _L_slices_size(L,D,l) = binomial(D-L+l,  D-L+1)
"""
_L_slices_size(L,D,l) = binomial(D-L+l,  D-L+1)

"""
    _sub_multi_indices!(sub_ids, α)

Given a positive multi-index `α`, sets in place in `sub_ids` the couples
(`id`, `d`) with `d` in 1:`N` for which the multi-index `αd⁻` = `α`-e`d` is
positive (that is `α`[`d`]>0), and `id` is the linear index of `αd⁻`
(see [`bernstein_term_id`](@ref)).

The function returns the number of sub indices set.
"""
function _sub_multi_indices!(sub_ids, α, ::Val{N}) where N
  @check length(sub_ids) >= N
  nb_sα = 0
  for i in 1:N
    α⁻ =  ntuple(k -> α[k]-Int(k==i), Val(N))
    if all(α⁻ .≥ 0)
      nb_sα += 1
      id⁻ = bernstein_term_id(α⁻)
      sub_ids[nb_sα] = (id⁻, i)
    end
  end
  return nb_sα
end

"""
    _sub_sub_multi_indices!(sub_ids, α, ::Val{N})

Like [`_sub_multi_indices`](@ref), but sets the triples (`id`, `t`, `q`) in `sub_ids`,
with `t,q` in 1:`N` for which the multi-index `αd⁻⁻` = `α`-e`t`-e`q` is positive,
and returns the number of triples set.
"""
function _sub_sub_multi_indices!(sub_ids, α, ::Val{N}) where N
  @check length(sub_ids) >= binomial(N+1,2)
  nb_ssα = 0
  for i in 1:N
    for j in i:N
      α⁻⁻ =  ntuple(k -> α[k]-Int(k==i)-Int(k==j), Val(N))
      if all(α⁻⁻ .≥ 0)
        nb_ssα += 1
        id⁻⁻ = bernstein_term_id(α⁻⁻)
        sub_ids[nb_ssα] = (id⁻⁻, i, j)
      end
    end
  end
  return nb_ssα
end

"""
    _sup_multi_indices!(sup_ids, α, ::Val{N})

Like [`_sub_multi_indices!`](@ref), but sets the indices for the `N` multi-indices
`αd⁺` = `α`+e`d` for 1≤d≤`N`, and returns `N`
"""
function _sup_multi_indices!(sup_ids, α, ::Val{N}) where N
  @check length(sup_ids) >= N
  for i in 1:N
    α⁺ =  ntuple(k -> α[k]+Int(k==i), Val(N))
    id⁺ = bernstein_term_id(α⁺)
    sup_ids[i] = (id⁺, i)
  end
  return N
end


################################
# nD evaluation implementation #
################################

# Overload _return_cache and _setsize for in place D-dimensional de Casteljau algorithm
function _return_cache(
  b::BernsteinBasisOnSimplex{D}, x,::Type{G},::Val{N_deriv}) where {D,G,N_deriv}

  T = eltype(G)
  K = get_order(b)
  np = length(x)
  ndof = length(b)
  ndof_scalar = binomial(K+D,D)

  r = CachedArray(zeros(G,(np,ndof)))
  s = MArray{Tuple{Vararg{D,N_deriv}},T}(undef)
  c = CachedVector(zeros(T,ndof_scalar))
  # The cache c here holds all scalar nD-Bernstein polynomials, no other caches needed for derivatives
  t = ntuple( _ -> nothing, Val(N_deriv))
  (r, s, c, t...)
end

function _setsize!(b::BernsteinBasisOnSimplex{D}, np, r, t...) where D
  K = get_order(b)
  ndof = length(b)
  ndof_scalar = binomial(K+D,D)
  setsize!(r,(np,ndof))
  setsize!(t[1],(ndof_scalar,))
end

function _evaluate_nd!(
  b::BernsteinBasisOnSimplex{D,V}, x,
  r::AbstractMatrix, i,
  c::AbstractVector{T}, VK::Val) where {D,V,T}

  λ = _cart_to_bary(x, b.cart_to_bary_matrix)
  c[1] = one(T)
  _downwards_de_Casteljau_nD!(c,λ,VK,Val(D))

  k = 1
  for s in c
    k = _cartprod_set_value!(r,i,s,k)
  end
end

function _gradient_nd!(
  b::BernsteinBasisOnSimplex{D,V}, x,
  r::AbstractMatrix{G}, i,
  c::AbstractVector{T},
  g::Nothing,
  s::MVector{D,T}, ::Val{K}) where {D,V,G,T,K}

  x_to_λ = b.cart_to_bary_matrix
  λ = _cart_to_bary(x, x_to_λ)

  c[1] = one(T)
  _downwards_de_Casteljau_nD!(c,λ,Val(K-1),Val(D))

  _grad_Bα_from_Bαm!(r,i,c,s,Val(K),Val(D),V,x_to_λ)
end

function _hessian_nd!(
  b::BernsteinBasisOnSimplex{D,V}, x,
  r::AbstractMatrix{G}, i,
  c::AbstractVector{T},
  g::Nothing,
  h::Nothing,
  s::MMatrix{D,D,T}, ::Val{K}) where {D,V,G,T,K}

  x_to_λ = b.cart_to_bary_matrix
  λ = _cart_to_bary(x, x_to_λ)

  c[1] = one(T)
  _downwards_de_Casteljau_nD!(c,λ,Val(K-2),Val(D))

  _hess_Bα_from_Bαmm!(r,i,c,s,Val(K),Val(D),V,x_to_λ)
end

# @generated functions as otherwise the time for computing the indices are the bottlneck...
@doc """
    _downwards_de_Casteljau_nD!(c, λ,::Val{K},::Val{D},::Val{K0}=Val(1))

Iteratively applies de Casteljau algorithm in reverse in place using `λ`s as
coefficients.

If `K0 = 1`, `λ` are the barycentric coordinates of some point `x` and `c[1] = 1`,
this computes all order `K` basis Bernstein polynomials at `x`:

`c[α_id] = B_α(x)  ∀α ∈ bernstein_terms(K,D)`

where `α_id` = [`bernstein_term_id`](@ref)(α).
"""
@generated function _downwards_de_Casteljau_nD!(c, λ,::Val{K},::Val{D},::Val{K0}=Val(1)) where {K,D,K0}
  z = zero(eltype(c))
  ex_v = Vector{Expr}()
  sub_ids = MVector{D+1,Tuple{Int,Int}}(undef)

  for Ki in K0:K
    # Iterations are  in reverse lexicographic order (left to right), because α-ei is
    # always stored on the left of α (as α-ei < α in lexicographic order), so the
    # erased B_β replaced by B_α won't be used to compute the remainings B_γ for |γ|=`K`
    # with γ>α in lexicographic order.
    terms = bernstein_terms(Ki,D)
    for (id,α) in Iterators.reverse(enumerate(terms)) # For all |α| = Ki
      # s = 0.
      push!(ex_v, :(s = $z))

      # For all |β| = |α|-1; β ≥ 0
      nb_sα = _sub_multi_indices!(sub_ids, α, Val(D+1))
      for (id_β, d) in take(sub_ids, nb_sα)
        # s +=  λ_d * B_β
        push!(ex_v, :(@inbounds s += λ[$d]*c[$id_β]))
      end

      # c[id] = B_α
      push!(ex_v, :(@inbounds c[$id] = s))
    end
  end
  return Expr(:block, ex_v...)
end

"""
    _de_Casteljau_nD!(c, λ,::Val{K},::Val{D},::Val{Kf}=Val(0))

Iteratively applies de Casteljau algorithm in place using `λ`s as
coefficients.

If `Kf = 0`, `λ` are the barycentric coordinates of some point `x` and `c` contains
the Bernstein coefficients ``c_α`` of a polynomial ``p`` (that is ``p(x) = ∑_α c_α B_α(x)`` for
``α`` in [`bernstein_terms`](@ref)(`K`,`D`) ), this computes

``c[1] = p(x)``

where the ``c_α`` must be initially stored in `c`[`α_id`], where
`α_id` = [`bernstein_term_id`](@ref)(α).
"""
function _de_Casteljau_nD! end #

@generated function _de_Casteljau_nD!(c, λ,::Val{K},::Val{D},::Val{Kf}=Val(0)) where {K,D,Kf}
  z = zero(eltype(c))
  ex_v = Vector{Expr}()
  sup_ids = MVector{D+1,Tuple{Int,Int}}(undef)

  for Ki in (K-1):-1:Kf
    # Iterations are in lexicographic order (right to left), because α+ei is
    # always stored on the right of α (as α+ei > α in lexicographic order), so the
    # erased B_β replaced by B_α won't be used to compute the remainings B_γ for |γ|=`K`
    # with γ<α in lexicographic order.
    terms = bernstein_terms(Ki,D)
    for (id,α) in enumerate(terms) # For all |α| = Ki
      # s = 0.
      push!(ex_v, :(s = $z))

      # For all |β| = |α|+1
      _sup_multi_indices!(sup_ids, α, Val(D+1))
      for (id_β, d) in sup_ids
        # s += λ_d * B_β
        push!(ex_v, :(@inbounds s += λ[$d]*c[$id_β]))
      end

      # c[id] = B_α (= s)
      push!(ex_v, :(@inbounds c[$id] = s))
    end
  end
  return Expr(:block, ex_v...)
end


# ∂q(B_α) = K ∑_{1 ≤ i ≤ N} ∂q(λi) B_β
# for  1 ≤ q ≤ D and β = α-ei
@generated function _grad_Bα_from_Bαm!(
    r,i,c,s,::Val{K},::Val{D},::Type{V},x_to_λ=nothing) where {K,D,V}

  ex_v = Vector{Expr}()
  ncomp = num_indep_components(V)
  z = zero(eltype(c))
  δ(i,j) = Int(i==j)
  sub_ids = MVector{D+1,Tuple{Int,Int}}(undef)

  for (id_α,α) in enumerate(bernstein_terms(K,D))
    push!(ex_v, :(@inbounds s .= $z))  # s = 0

    nb_sα = _sub_multi_indices!(sub_ids, α, Val(D+1))
    for (id_β, i) in take(sub_ids, nb_sα) # β = α - ei
      push!(ex_v, :(@inbounds B_β = c[$id_β]))
      # s[q] = Σ_β ∂q(λi) B_β
      for q in 1:D
        if x_to_λ == Nothing
          # ∇λ(eq)_i = δ_{q+1,i} - δ_1i
          Cqi = δ(i,q+1) - δ(1,i)
          iszero(Cqi) || push!(ex_v, :(@inbounds s[$q] += $Cqi*B_β))
        else
          # ∂q(λi) = ei (x_to_λ*(e1 - e{q+1}) - x_to_λ*(e1)) = ei * x_to_λ * e{q+1}
          # ∂q(λi) = x_to_λ[i,q+1]
          push!(ex_v, :(@inbounds s[$q] += x_to_λ[$i,$(q+1)]*B_β))
        end
      end
    end
    push!(ex_v, :(@inbounds s .*= $K)) # s = Ks.

    k = ncomp*(id_α-1) + 1
    push!(ex_v, :(_cartprod_set_derivative!(r,i,s,$k,V)))
  end

  return Expr(:block, ex_v...)
end

# ∂t(∂q(B_α)) = K(K-1) ∑_{1 ≤ i,j ≤ N} ∂t(λj) ∂q(λi) B_β
# for  1 ≤ t,q ≤ D and β = α - ei - ej
@generated function _hess_Bα_from_Bαmm!(
    r,i,c,s,::Val{K},::Val{D},::Type{V},x_to_λ=nothing) where {K,D,V}

  ex_v = Vector{Expr}()
  ncomp = num_indep_components(V)
  z = zero(eltype(c))
  δ(i,j) = Int(i==j)
  C(q,t,i,j) = (δ(i,q+1)-δ(i,1))*(δ(j,t+1)-δ(j,1))
  N_max_ssα = binomial(D+2,2)
  sub_sub_α_ids = MVector{N_max_ssα, NTuple{3,Int}}(undef)

  for (id_α,α) in enumerate(bernstein_terms(K,D))
    push!(ex_v, :(@inbounds s .= $z))     # s = 0

    nb_subsub = _sub_sub_multi_indices!(sub_sub_α_ids, α, Val(D+1))
    for (id_β, i, j) in take(sub_sub_α_ids, nb_subsub)

      push!(ex_v, :(@inbounds B_β = c[$id_β]))
      for t in 1:D
        for q in 1:D
          if x_to_λ == Nothing
          # s[t,q] = Σ_β B_β (δ_iq - δ_iN)*(δ_jt - δ_jN)
                   Cβ  = C(q,t,i,j)
            if i≠j Cβ += C(q,t,j,i) end
            iszero(Cβ) || push!(ex_v, :(@inbounds s[$t,$q] += $Cβ*B_β))
          else
                   push!(ex_v, :(@inbounds C =  x_to_λ[$i,$(q+1)]*x_to_λ[$j,$(t+1)]))
            if i≠j push!(ex_v, :(@inbounds C += x_to_λ[$j,$(q+1)]*x_to_λ[$i,$(t+1)])) end
            push!(ex_v, :(@inbounds s[$t,$q] += C*B_β))
          end
        end
      end
    end
    push!(ex_v, :(@inbounds s .*= $(K*(K-1))) )  # s = K(K-1)s

    k = ncomp*(id_α-1) + 1
    push!(ex_v, :(_cartprod_set_derivative!(r,i,s,$k,V)))
  end

  return Expr(:block, ex_v...)
end


#######################################################
#### Bernstein bases on simplices Naive implementation#
#######################################################
#
#  """
#      BernsteinBasisOnSimplex{D,V} <: PolynomialBasis{D,V,Bernstein}
#
#  This basis uses barycentric coordinates defined by the vertices of the
#  reference `D`-simplex.
#  """
#  struct BernsteinBasisOnSimplex{D,V} <: PolynomialBasis{D,V,Bernstein}
#    function BernsteinBasisOnSimplex{D}(::Type{V},order::Int) where {D,V}
#      K = Int(order)
#      new{D,V}()
#    end
#  end
#
#  function BernsteinBasisOnSimplex(::Val{D},::Type{V},order::Int) where {D,V}
#    BernsteinBasisOnSimplex{D}(V,order)
#  end
#
#  Base.size(::BernsteinBasisOnSimplex{D,V}) where {D,V} = (num_indep_components(V)*binomial(D+K,D),)
#  get_exponents(::BernsteinBasisOnSimplex{D,V}) where {D,V} = bernstein_terms(K,D)
#
#  ################################
#  # nD evaluation implementation #
#  ################################
#
#  # Overload _return_cache and _setsize to add +1 coordinate cache in t
#  function _return_cache(
#    f::BernsteinBasisOnSimplex{D}, x,::Type{G},::Val{N_deriv}) where {D,G,N_deriv}
#
#    T = eltype(G)
#    np = length(x)
#    ndof = length(f)
#    ndof_1d = get_order(f) + 1
#    r = CachedArray(zeros(G,(np,ndof)))
#    s = MArray{Tuple{Vararg{D,N_deriv}},T}(undef)
#    bernstein_D = D+1 # There are D+1 barycentric coordinates
#    t = ntuple( _ -> CachedArray(zeros(T,(bernstein_D,ndof_1d))), Val(N_deriv+1))
#    (r, s, t...)
#  end
#  function _setsize!(f::BernsteinBasisOnSimplex{D}, np, r, t...) where D
#    ndof = length(f)
#    ndof_1d = get_order(f) + 1
#    setsize!(r,(np,ndof))
#    bernstein_D = D+1 # There are D+1 barycentric coordinates
#    for c in t
#      setsize!(c,(bernstein_D,ndof_1d))
#    end
#  end
#
#
#  function _evaluate_nd!(
#    b::BernsteinBasisOnSimplex{D,V}, x,
#    r::AbstractMatrix{V}, i,
#    c::AbstractMatrix{T}) where {D,V,T}
#
#    K = get_order(b)
#    terms  = _get_terms(b)
#    coefs = multinoms(Val(K),Val(D))
#
#    λ = _cart_to_bary(x,nothing)
#
#    for d in 1:(D+1)
#      _evaluate_1d!(Monomial,Val(K),c,λ,d) # compute powers 0:K of all bary. coords.
#    end
#
#    k = 1
#    for (ci,m) in zip(terms,coefs)
#
#      for d in 1:(D+1)
#        @inbounds m *= c[d,ci[d]]
#      end
#
#      k = _cartprod_set_value!(r,i,m,k)
#    end
#  end
#
#  function _gradient_nd!(
#    b::BernsteinBasisOnSimplex{D,V}, x,
#    r::AbstractMatrix{G}, i,
#    c::AbstractMatrix{T},
#    g::AbstractMatrix{T},
#    s::MVector{D,T}) where {D,V,G,T}
#
#    K = get_order(b)
#    N = D+1
#    terms = _get_terms(b)
#    coefs = multinoms(Val(K),Val(D))
#
#    λ = _cart_to_bary(x,nothing)
#
#    for d in 1:N
#      _derivatives_1d!(Monomial,Val(K),(c,g),λ,d)
#    end
#
#    k = 1
#    @inbounds for (ci,m) in zip(terms,coefs)
#
#      for i in eachindex(s)
#        s[i] = m
#      end
#
#      for q in 1:D
#        for d in 1:D
#          if d != q
#            s[q] *= c[d,ci[d]]
#          else
#            s[q] *= g[q,ci[q]]*c[N,ci[N]] - g[N,ci[N]]*c[q,ci[q]]
#          end
#        end
#      end
#
#      k = _cartprod_set_derivative!(r,i,s,k,V)
#    end
#  end
#
#  function _hessian_nd!(
#    b::BernsteinBasisOnSimplex{D,V}, x,
#    r::AbstractMatrix{G}, i,
#    c::AbstractMatrix{T},
#    g::AbstractMatrix{T},
#    h::AbstractMatrix{T},
#    s::MMatrix{D,D,T}) where {D,V,G,T}
#
#    N = D+1
#    terms = _get_terms(b)
#    coefs = multinoms(Val(K),Val(D))
#
#    λ = _cart_to_bary(x,nothing)
#
#    for d in 1:N
#      _derivatives_1d!(Monomial,Val(K),(c,g,h),λ,d)
#    end
#
#    k = 1
#    @inbounds for (ci,m) in zip(terms,coefs)
#
#      for i in eachindex(s)
#        s[i] = m
#      end
#
#      for t in 1:D
#        for q in 1:D
#          for d in 1:D
#            if d != q && d != t
#              # if q == t, D-1 factors
#              # else,      D-2 factors
#              s[t,q] *= c[d,ci[d]]
#            elseif q == t # == d
#              # +2 factors -> D+1
#              s[t,q] *= (h[d,ci[d]]*c[N,ci[N]] -2g[d,ci[d]]*g[N,ci[N]] + c[d,ci[d]]*h[N,ci[N]])
#            elseif d == q # q ≠ t, we multiply once with the factors with q and t derivative terms
#              # +3 factors -> D+1
#              s[t,q] *=(  g[t,ci[t]]*g[q,ci[q]]*c[N,ci[N]]
#                        - g[t,ci[t]]*c[q,ci[q]]*g[N,ci[N]]
#                        - c[t,ci[t]]*g[q,ci[q]]*g[N,ci[N]]
#                        + c[t,ci[t]]*c[q,ci[q]]*h[N,ci[N]])
#            end
#          end
#        end
#      end
#
#      k = _cartprod_set_derivative!(r,i,s,k,V)
#    end
#  end
#
# """
#     multinoms(::Val{K}, ::Val{D})
#
# Returns the tuple of multinomial coefficients for each term in
# [`bernstein_terms`](@ref)(`K`,`D`). For e.g. a term `t`, the
# multinomial can be computed by `factorial(sum(t)) ÷ prod(factorial.(t)`
# """
# @generated function multinoms(::Val{K},::Val{D}) where {K,D}
#   terms = bernstein_terms(K,D)
#   multinomials = tuple( (multinomial(α...) for α in terms)... )
#   Meta.parse("return $multinomials")
# end
#
