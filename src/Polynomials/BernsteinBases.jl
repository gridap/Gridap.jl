"""
    Bernstein <: Polynomial

Type representing Bernstein polynomials, c.f. [Bernstein polynomials](@ref) section.
"""
struct Bernstein <: Polynomial end

isHierarchical(::Type{Bernstein}) = false


##########################################
# Uniform tensor product Bernstein bases #
##########################################

"""
    BernsteinBasis{D,V,K} = UniformPolyBasis{D,V,K,Bernstein}

Alias for uniform Bernstein multivariate scalar' or `Multivalue`'d basis.
"""
const BernsteinBasis{D,V,K} = UniformPolyBasis{D,V,K,Bernstein}

"""
    BernsteinBasis(::Val{D}, ::Type{V}, order::Int, terms::Vector)
    BernsteinBasis(::Val{D}, ::Type{V}, order::Int [, filter::Function])
    BernsteinBasis(::Val{D}, ::Type{V}, orders::Tuple [, filter::Function])

High level constructors of [`BernsteinBasis`](@ref).
"""
BernsteinBasis(args...) = UniformPolyBasis(Bernstein, args...)


################################
# Bernstein bases on simplices #
################################

"""
    BernsteinBasisOnSimplex{D,V,K} <: PolynomialBasis{D,V,K,Bernstein}

This basis uses barycentric coordinates defined by the vertices of the
reference `D`-simplex.
"""
struct BernsteinBasisOnSimplex{D,V,K} <: PolynomialBasis{D,V,K,Bernstein}
  function BernsteinBasisOnSimplex{D}(::Type{V},order::Int) where {D,V}
    K = Int(order)
    new{D,V,K}()
  end
end

function BernsteinBasisOnSimplex(::Val{D},::Type{V},order::Int) where {D,V}
  BernsteinBasisOnSimplex{D}(V,order)
end

Base.size(::BernsteinBasisOnSimplex{D,V,K}) where {D,V,K} = (num_indep_components(V)*binomial(D+K,D),)
_get_terms(::BernsteinBasisOnSimplex{D,V,K}) where {D,V,K} = bernstein_terms(Val(K), Val(D))

function get_exponents(b::BernsteinBasisOnSimplex)
  indexbase = 1
  terms = _get_terms(b)
  [Tuple(t) .- indexbase for t in terms]
end


###########
# Helpers #
###########

"""
    _cart_to_bary(x::Point{D,T})

Converts the cartesian coordinates `x` into the barycentric coordinates with
respect to the reference simplex, that is `λ`=(x1, ..., xD, 1-x1-x2-...-xD).
"""
@inline function _cart_to_bary(x::Point{D,T}) where {D,T}
  λ = zero(MVector{D+1,T})

  s = zero(T)
  @inbounds begin
    for d in 1:D
      xd = x[d]
      s += xd
      λ[d] = x[d]
    end
    λ[D+1] = 1-s
  end

  return Tuple(λ)
end

"""
    bernstein_terms(::Val{K},::Val{D})

Return the set of multi-indices for the `D`-dimensional Bernstein basis of
order `K`, that is:

TODO
"""
@generated function bernstein_terms(::Val{K},::Val{D}) where {K,D}
  indexbase = 1
  multi_exponents = collect( tuple((v .+ indexbase)...) for v in multiexponents(D+1,K))
  terms = tuple(multi_exponents...)
  Meta.parse("return $terms")
end

"""
    binom(::Val{K}, ::Val{I})

Returns the binomial coefficient C(K,I).
"""
_binomial(::Val{K},::Val{I}) where {K,I} = binomial(K,I)

"""
    binoms(::Val{K})

Returns the tuple of binomials ( C₍ₖ₀₎, C₍ₖ₁₎, ..., C₍ₖₖ₎ ).
"""
binoms(::Val{K}) where K = ntuple( i -> binomial(K,i-1), Val(K+1))

"""
    multinoms(::Val{K}, ::Val{D})

Returns the tuple of multinomial coefficients for each term in
[`bernstein_terms`](@ref)(Val(`K`),Val(`D`)). For e.g. a term `t`, the
multinomial can be computed by `factorial(sum(t)) ÷ prod(factorial.(t)`
"""
@generated function multinoms(::Val{K},::Val{D}) where {K,D}
  terms = bernstein_terms(Val(K),Val(D))
  indexbase = 1
  multinomials = tuple( (multinomial((t .- indexbase)...) for t in terms)... )
  Meta.parse("return $multinomials")
end


################################
# nD evaluation implementation #
################################

# Overlead _return_cache and _setsize to add +1 coordinate cache in t
function _return_cache(
  f::BernsteinBasisOnSimplex{D}, x,::Type{G},::Val{N_deriv}) where {D,G,N_deriv}

  @assert D == length(eltype(x)) "Incorrect number of point components"
  T = eltype(G)
  np = length(x)
  ndof = length(f)
  ndof_1d = get_order(f) + 1
  r = CachedArray(zeros(G,(np,ndof)))
  bernstein_D = D+1 # There are D+1 barycentryc coordinates
  t = ntuple( _ -> CachedArray(zeros(T,(bernstein_D,ndof_1d))), Val(N_deriv+1))
  (r, t...)
end
function _setsize!(f::BernsteinBasisOnSimplex{D}, np, r, t...) where D
  ndof = length(f)
  ndof_1d = get_order(f) + 1
  setsize!(r,(np,ndof))
  bernstein_D = D+1 # There are D+1 barycentryc coordinates
  for c in t
    setsize!(c,(bernstein_D,ndof_1d))
  end
end


function _evaluate_nd!(
  b::BernsteinBasisOnSimplex{D,V,K}, x,
  r::AbstractMatrix{V}, i,
  c::AbstractMatrix{T}) where {D,V,K,T}

  terms  = _get_terms(b)
  coefs = multinoms(Val(K),Val(D))

  λ = _cart_to_bary(x)

  for d in 1:(D+1)
    _evaluate_1d!(Monomial,Val(K),c,λ,d) # compute powers 0:K of all bary. coords.
  end

  k = 1
  for (ci,m) in zip(terms,coefs)

    for d in 1:(D+1)
      @inbounds m *= c[d,ci[d]]
    end

    k = _uniform_set_value!(r,i,m,k)
  end
end

function _gradient_nd!(
  b::BernsteinBasisOnSimplex{D,V,K}, x,
  r::AbstractMatrix{G}, i,
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  s::MVector{D,T}) where {D,V,K,G,T}

  N = D+1
  terms = _get_terms(b)
  coefs = multinoms(Val(K),Val(D))

  λ = _cart_to_bary(x)

  for d in 1:N
    _derivatives_1d!(Monomial,Val(K),(c,g),λ,d)
  end

  k = 1
  @inbounds for (ci,m) in zip(terms,coefs)

    for i in eachindex(s)
      s[i] = m
    end

    for q in 1:D
      for d in 1:D
        if d != q
          s[q] *= c[d,ci[d]]
        else
          s[q] *= g[q,ci[q]]*c[N,ci[N]] - g[N,ci[N]]*c[q,ci[q]]
        end
      end
    end

    k = _uniform_set_derivative!(r,i,s,k,V)
  end
end

function _hessian_nd!(
  b::BernsteinBasisOnSimplex{D,V,K}, x,
  r::AbstractMatrix{G}, i,
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  h::AbstractMatrix{T},
  s::MMatrix{D,D,T}) where {D,V,K,G,T}

  N = D+1
  terms = _get_terms(b)
  coefs = multinoms(Val(K),Val(D))

  λ = _cart_to_bary(x)

  for d in 1:N
    _derivatives_1d!(Monomial,Val(K),(c,g,h),λ,d)
  end

  k = 1
  @inbounds for (ci,m) in zip(terms,coefs)

    for i in eachindex(s)
      s[i] = m
    end

    for t in 1:D
      for q in 1:D
        for d in 1:D
          if d != q && d != t
            # if q == t, D-1 factors
            # else,      D-2 factors
            s[t,q] *= c[d,ci[d]]
          elseif q == t # == d
            # +2 factors -> D+1
            s[t,q] *= (h[d,ci[d]]*c[N,ci[N]] -2g[d,ci[d]]*g[N,ci[N]] + c[d,ci[d]]*h[N,ci[N]])
          elseif d == q # q ≠ t, we multiply once with the factors with q and t derivative terms
            # +3 factors -> D+1
            s[t,q] *=(  g[t,ci[t]]*g[q,ci[q]]*c[N,ci[N]]
                      - g[t,ci[t]]*c[q,ci[q]]*g[N,ci[N]]
                      - c[t,ci[t]]*g[q,ci[q]]*g[N,ci[N]]
                      + c[t,ci[t]]*c[q,ci[q]]*h[N,ci[N]])
          end
        end
      end
    end

    k = _uniform_set_derivative!(r,i,s,k,V)
  end
end


################################
# 1D evaluation implementation #
################################

function _evaluate_1d!(::Type{Bernstein},::Val{0},v::AbstractMatrix{T},x,d) where {T<:Number}
  @inbounds v[d,1] = one(T)
end

@inline function _De_Casteljau_step_1D!(v,d,i,λ1,λ2)
  # i = k+1

  # vₖ <- xvₖ₋₁            # Bᵏₖ(x) = x*Bᵏ⁻¹ₖ₋₁(x)
  v[d,i] = λ2*v[d,i-1]
  # vⱼ <- xvⱼ₋₁ + (1-x)vⱼ  # Bᵏⱼ(x) = x*Bᵏ⁻¹ⱼ₋₁(x) + (1-x)*Bᵏ⁻¹ⱼ(x) for j = k-1, k-2, ..., 1
  for l in i-1:-1:2
    v[d,l] = λ2*v[d,l-1] + λ1*v[d,l]
  end
  # v₀ <- (1-x)v₀          # Bᵏ₀(x) = (1-x)*Bᵏ⁻¹₀(x)
  v[d,1] = λ1*v[d,1]
end

# jth Bernstein poly of order K at x:
# Bᵏⱼ(x) = binom(K,j) * x^j * (1-x)^(K-j) = x*Bᵏ⁻¹ⱼ₋₁(x) + (1-x)*Bᵏ⁻¹ⱼ(x)
function _evaluate_1d!(::Type{Bernstein},::Val{K},v::AbstractMatrix{T},x,d) where {K,T<:Number}
  @inbounds begin
    n = K + 1 # n > 1
    λ2 = x[d]
    λ1 = one(T) - λ2

    # In place De Casteljau: init with B¹₀(x)=x and B¹₁(x)=1-x
    v[d,1] = λ1
    v[d,2] = λ2

    for i in 3:n
      _De_Casteljau_step_1D!(v,d,i,λ1,λ2)
    end
  end
  # still optimisable for K > 2/3:
  # - compute bj = binoms(k,j) at compile time (binoms(Val(K)) function)
  # - compute vj = xʲ*(1-x)ᴷ⁻ʲ recursively in place like De Casteljau (saving half the redundant multiplications)
  # - do it in a stack allocated cache (MVector, Bumber.jl)
  # - @simd affect bj * vj in v[d,i] for all j
end

function _gradient_1d!(::Type{Bernstein},::Val{0},g::AbstractMatrix{T},x,d) where {T<:Number}
  @inbounds g[d,1] = zero(T)
end
function _gradient_1d!(::Type{Bernstein},::Val{1},g::AbstractMatrix{T},x,d) where {T<:Number}
  o = one(T)
  @inbounds g[d,1] = -o
  @inbounds g[d,2] =  o
end

# First derivative of the jth Bernstein poly of order K at x:
# (Bᵏⱼ)'(x) = K * ( Bᵏ⁻¹ⱼ₋₁(x) - Bᵏ⁻¹ⱼ(x) )
#           = K * x^(j-1) * (1-x)^(K-j-1) * ((1-x)*binom(K-1,j-1) - x*binom(K-1,j))
function _gradient_1d!(::Type{Bernstein},::Val{K}, g::AbstractMatrix{T},x,d) where {K,T<:Number}
  @inbounds begin
    n = K + 1 # n > 2

    # De Casteljau for Bᵏ⁻¹ⱼ for j = k-1, k-2, ..., 1
    _evaluate_1d!(Bernstein,Val(K-1),g,x,d)

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


function _hessian_1d!(::Type{Bernstein},::Val{0},h::AbstractMatrix{T},x,d) where {T<:Number}
  @inbounds h[d,1] = zero(T)
end
function _hessian_1d!(::Type{Bernstein},::Val{1},h::AbstractMatrix{T},x,d) where {T<:Number}
  @inbounds h[d,1] = zero(T)
  @inbounds h[d,2] = zero(T)
end
function _hessian_1d!(::Type{Bernstein},::Val{2},h::AbstractMatrix{T},x,d) where {T<:Number}
  o = one(T)
  @inbounds h[d,1] =  2o
  @inbounds h[d,2] = -4o
  @inbounds h[d,3] =  2o
end

# Second derivative of the jth Bernstein poly of order K at x:
# (Bᵏⱼ)''(x) = K(K-1) * ( Bᵏ⁻²ⱼ₋₂(x) -2*Bᵏ⁻²ⱼ₋₁(x) + Bᵏ⁻²ⱼ(x) )
#            = K(K-1) * x^(j-2) * (1-x)^(K-j-2) * ( (1-x)^2*binom(K-2,j-2)
#                  - 2x*(1-x)*binom(K-2,j-1) + (x)^2*binom(K-2,j)
#              )
function _hessian_1d!(::Type{Bernstein},::Val{K},h::AbstractMatrix{T},x,d) where {K,T<:Number}
  @inbounds begin
    n = K + 1 # n > 3
    KK = K*(K-1)

    # De Casteljau for Bᵏ⁻²ⱼ for j = k-2, k-3, ..., 1
    _evaluate_1d!(Bernstein,Val(K-2),h,x,d)

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

function _derivatives_1d!(::Type{Bernstein},v::Val_01,t::NTuple{2},x,d)
  @inline _evaluate_1d!(Bernstein, v, t[1], x, d)
  @inline _gradient_1d!(Bernstein, v, t[2], x, d)
end

function _derivatives_1d!(::Type{Bernstein},::Val{K},t::NTuple{2},x,d) where K
  @inbounds begin
    n = K + 1 # n > 2
    v, g = t

    λ2 = x[d]
    λ1 = one(eltype(v)) - λ2

    # De Casteljau for Bᵏ⁻¹ⱼ for j = k-1, k-2, ..., 1
    _evaluate_1d!(Bernstein,Val(K-1),v,x,d)

    # Compute gradients as _gradient_1d!
    g[d,n] = K*v[d,n-1]
    @simd for l in n-1:-1:2
      g[d,l] = K*(v[d,l-1] - v[d,l])
    end
    g[d,1] = -K*v[d,1]

    # Last step of De Casteljau for _evaluate_1d!
    _De_Casteljau_step_1D!(v,d,n,λ1,λ2)
  end
end

function _derivatives_1d!(::Type{Bernstein},v::Val_012,t::NTuple{3},x,d)
  @inline _evaluate_1d!(Bernstein, v, t[1], x, d)
  @inline _gradient_1d!(Bernstein, v, t[2], x, d)
  @inline _hessian_1d!( Bernstein, v, t[3], x, d)
end

function _derivatives_1d!(::Type{Bernstein},::Val{K},t::NTuple{3},x,d) where K
  @inbounds begin
    n = K + 1 # n > 3
    v, g, h = t

    λ2 = x[d]
    λ1 = one(eltype(v)) - λ2

    # De Casteljau until Bᵏ⁻²ⱼ ∀j
    _evaluate_1d!(Bernstein,Val(K-2),v,x,d)

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
    _De_Casteljau_step_1D!(v,d,n-1,λ1,λ2)

    # Compute gradients as in _gradient_1d!
    g[d,n] = K*v[d,n-1]
    @simd for l in n-1:-1:2
      g[d,l] = K*(v[d,l-1] - v[d,l])
    end
    g[d,1] = -K*v[d,1]

    # Last step of De Casteljau for _evaluate_1d!
    _De_Casteljau_step_1D!(v,d,n,λ1,λ2)
  end
end


###################################################
# Bernstein bases on simplices using de Casteljau #
###################################################

"""
    BernsteinBasisOnSimplexDC{D,V,K} <: PolynomialBasis{D,V,K,Bernstein}

This basis uses barycentric coordinates defined by the vertices of the
reference `D`-simplexDC.
"""
struct BernsteinBasisOnSimplexDC{D,V,K} <: PolynomialBasis{D,V,K,Bernstein}
  function BernsteinBasisOnSimplexDC{D}(::Type{V},order::Int) where {D,V}
    K = Int(order)
    new{D,V,K}()
  end
end

function BernsteinBasisOnSimplexDC(::Val{D},::Type{V},order::Int) where {D,V}
  BernsteinBasisOnSimplexDC{D}(V,order)
end

Base.size(::BernsteinBasisOnSimplexDC{D,V,K}) where {D,V,K} = (num_indep_components(V)*binomial(D+K,D),)
_get_terms(::BernsteinBasisOnSimplexDC{D,V,K}) where {D,V,K} = bernstein_terms(Val(K), Val(D))

function get_exponents(b::BernsteinBasisOnSimplexDC)
  indexbase = 1
  terms = _get_terms(b)
  [Tuple(t) .- indexbase for t in terms]
end


########################
# de Casteljau helpers #
########################

"""
    _L_slice(L,α::NTuple{N}) where N = sum(last(α,N-L))

where `L` ∈ 1:N

For a given positive Bernstein term `α`, return the index (starting
from 0) of the (D-`L`)-slice to which `α` belongs within the (D-`L`-1)-slice of
the D-multiexponent simplex (D = `N`-1).

In a D-multiexponent simplex of elements `α`, flattened in a vector in
lexicographic order, the (D-`L`)-slices are the consecutive `α`s
having iddentical first `L` indices `α`.

For example, the (3-1)=2 slices of the tetrahedral multiexponents (3 simplex) are triangle multiexponents.
"""
_L_slice(L,α::NTuple{N}) where N = sum(last(α,N-L))

"""
    _L_slices_size(L,D,l) = binomial(D-L+l,  D-L+1)

Return the length of the `l`-1 first (`D`-`L`)-slices in the
`D`-multiexponent simplex (flattened in lexicographic order).
Those numbers are the "(`D`-`L`)-simplex numbers".
"""
_L_slices_size(L,D,l) = binomial(D-L+l,  D-L+1)

"""
    _simplex_multi_id_to_linear_id(α::NTuple{N,Int})

For a given positive Bernstein term `α`, return the linear index of
the associated Bernstein polynomial, that is the `i` such that

    (i,α) ∈ enumerate(bernstein_terms(Val(K),Val(N-1))

where K = sum(`α`).
"""
function _simplex_multi_id_to_linear_id(α::NTuple{N}) where N
  D = N-1
  i = sum( _L_slices_size(L, D, _L_slice(L,α)) for L in 1:D) + 1
  i
end

"""
    _sub_multi_indices(α::NTuple{N,Int})

Given a positive multi-index `α`, return a tuple of couples
(`id`, `d`) for `d` in 1:`N` for which the multi-index `αd⁻` = `α`-e`d` is
positive (that is `α`[`d`]>0), and `id` is the linear index of `αd⁻`
(see [`_simplex_multi_id_to_linear_id`](@ref)).
"""
function _sub_multi_indices(α::NTuple{N,Int}) where N
  sub_ids = tuple()
  for d in 1:N
    if α[d] > 0
      αd⁻ = ntuple( j -> α[j] - ==(d,j), N)
      id⁻ = _simplex_multi_id_to_linear_id(αd⁻)
      sub_ids = (sub_ids..., (id⁻, d, αd⁻))
    end
  end
  return sub_ids
end

function _de_Casteljau_indices(K,D)
  terms = bernstein_terms(Val(K),Val(D))
  # The Iterators.reverse is important to avoid erasing values of previous
  # iteration that are still needed later (for the in-place de Casteljau iteration).
  rev_terms = Iterators.reverse(enumerate(terms))
  return ( (id,_sub_multi_indices(α .- 1)) for (id,α) in rev_terms )
end

################################
# nD evaluation implementation #
################################

# Overlead _return_cache and _setsize for in place D-dimensional de Casteljau algorithm
function _return_cache(
  #f::BernsteinBasisOnSimplexDC{D}, x,::Type{G},::Val{N_deriv}) where {D,G,N_deriv}
  f::BernsteinBasisOnSimplexDC{D}, x,::Type{G},::Val{0}) where {D,G}
  N_deriv = 0

  @assert D == length(eltype(x)) "Incorrect number of point components"
  T = eltype(G)
  K = get_order(f)
  np = length(x)
  ndof = length(f)
  ndof_scalar = _binomial(Val(K+D),Val(D))

  r = CachedArray(zeros(G,(np,ndof)))
  c = CachedVector(zeros(T,ndof_scalar))
  # The cache t here holds all scalar nD-Bernstein polymials, no other caches needed for derivatives
  t = ntuple( i -> isone(i) ? c : nothing , Val(N_deriv+1))
  (r, t...)
end
#function _setsize!(f::BernsteinBasisOnSimplexDC{D}, np, r, t...) where D
function _setsize!(f::BernsteinBasisOnSimplexDC{D}, np, r, t) where D
  K = get_order(f)
  ndof = length(f)
  ndof_scalar = _binomial(Val(K+D),Val(D))
  setsize!(r,(np,ndof))
  #setsize!(t[1],(ndof_scalar,))
  setsize!(t,(ndof_scalar,))
end

const Val_12 = Union{Val{1},Val{2}}
function _return_cache(
  f::BernsteinBasisOnSimplexDC{D}, x,::Type{G},val::Val_12) where {D,G}

  @assert D == length(eltype(x)) "Incorrect number of point components"
  T = eltype(G)
  np = length(x)
  ndof = length(f)
  ndof_1d = get_order(f) + 1
  r = CachedArray(zeros(G,(np,ndof)))
  bernstein_D = D+1 # There are D+1 barycentryc coordinates
  t = ntuple( _ -> CachedArray(zeros(T,(bernstein_D,ndof_1d))), Val(typeof(val).parameters[1]+1))
  (r, t...)
end
function _setsize!(f::BernsteinBasisOnSimplexDC{D}, np, r, t...) where D
  ndof = length(f)
  ndof_1d = get_order(f) + 1
  setsize!(r,(np,ndof))
  bernstein_D = D+1 # There are D+1 barycentryc coordinates
  for c in t
    setsize!(c,(bernstein_D,ndof_1d))
  end
end

# @generated function as otherwise the time and allocation for
# computing the indices are the bottlneck...
@generated function _De_Casteljau_nD!(c, λ, ::Val{K}, ::Val{D}) where {K,D}
  z = zero(eltype(c))
  ex_v = Vector{Expr}()
  for Ki in 1:K
    for (id,sub_ids) in _de_Casteljau_indices(Ki,D)

      # s = 0.
      push!(ex_v, :(s = $z))
      for (idα⁻, d) in sub_ids
        # s += λ[d]*c[idα⁻]
        push!(ex_v, :(@inbounds s += λ[$d]*c[$idα⁻]))
      end

      # c[id] = s
      push!(ex_v, :(@inbounds c[$id] = s))
    end
  end
  return Expr(:block, ex_v...)
end

function _evaluate_nd!(
  b::BernsteinBasisOnSimplexDC{D,V,K}, x,
  r::AbstractMatrix{V}, i,
  c::AbstractVector{T}) where {D,V,K,T}

  λ = _cart_to_bary(x)
  c[1] = one(T)
  _De_Casteljau_nD!(c,λ,Val(K),Val(D))

  k = 1
  for s in c
    k = _uniform_set_value!(r,i,s,k)
  end
end

function _gradient_nd!(
  b::BernsteinBasisOnSimplexDC{D,V,K}, x,
  r::AbstractMatrix{G}, i,
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  s::MVector{D,T}) where {D,V,K,G,T}

  N = D+1
  terms = _get_terms(b)
  coefs = multinoms(Val(K),Val(D))

  λ = _cart_to_bary(x)

  for d in 1:N
    _derivatives_1d!(Monomial,Val(K),(c,g),λ,d)
  end

  k = 1
  @inbounds for (ci,m) in zip(terms,coefs)

    for i in eachindex(s)
      s[i] = m
    end

    for q in 1:D
      for d in 1:D
        if d != q
          s[q] *= c[d,ci[d]]
        else
          s[q] *= g[q,ci[q]]*c[N,ci[N]] - g[N,ci[N]]*c[q,ci[q]]
        end
      end
    end

    k = _uniform_set_derivative!(r,i,s,k,V)
  end
end

function _hessian_nd!(
  b::BernsteinBasisOnSimplexDC{D,V,K}, x,
  r::AbstractMatrix{G}, i,
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  h::AbstractMatrix{T},
  s::MMatrix{D,D,T}) where {D,V,K,G,T}

  N = D+1
  terms = _get_terms(b)
  coefs = multinoms(Val(K),Val(D))

  λ = _cart_to_bary(x)

  for d in 1:N
    _derivatives_1d!(Monomial,Val(K),(c,g,h),λ,d)
  end

  k = 1
  @inbounds for (ci,m) in zip(terms,coefs)

    for i in eachindex(s)
      s[i] = m
    end

    for t in 1:D
      for q in 1:D
        for d in 1:D
          if d != q && d != t
            # si q == t, D-1 facteurs
            # sinon,     D-2 facteurs
            s[t,q] *= c[d,ci[d]]
          elseif q == t # == d
            # +2 facteurs -> D+1
            s[t,q] *= (h[d,ci[d]]*c[N,ci[N]] -2g[d,ci[d]]*g[N,ci[N]] + c[d,ci[d]]*h[N,ci[N]])
          elseif d == q # q ≠ t, we multiply once with the factors with q and t derivative terms
            # +3 facteurs -> D+1
            s[t,q] *=(  g[t,ci[t]]*g[q,ci[q]]*c[N,ci[N]]
                      - g[t,ci[t]]*c[q,ci[q]]*g[N,ci[N]]
                      - c[t,ci[t]]*g[q,ci[q]]*g[N,ci[N]]
                      + c[t,ci[t]]*c[q,ci[q]]*h[N,ci[N]])
          end
        end
      end
    end

    k = _uniform_set_derivative!(r,i,s,k,V)
  end
end


