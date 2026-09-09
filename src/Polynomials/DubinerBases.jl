############################################################################################
# The Dubiner / Koornwinder basis: the L²-orthonormal polynomial basis of `Pᴷ` on
# the reference simplex, in any dimension.
#
# WHY
#
# `LegendreBasis` is a *tensor product* of 1D Legendre polynomials, so it is
# L²-orthogonal on a segment and on an n-cube — and on nothing else. On the
# reference triangle its non-constant members have nonzero mean (-1/√3 for the
# linears), so "degree > q" and "orthogonal to P_q" do not coincide there. With
# an orthogonal basis on the simplex they do, and a complement can be selected by
# a filter:
#
#     DubinerBasis(Val(2), T, 2, _pk_minus_pq_filter(2, 0))   # P₂(K) ∩ P₀(K)^⊥
#
# exactly as edge weights already do with `LegendreBasis` on a segment. That is
# what the augmented (constrained-prebasis) reference FEs need for their cell
# moment weights.
#
# THE BASIS
#
# On the reference `D`-simplex `{x : xᵢ ≥ 0, Σxᵢ ≤ 1}`, for a multi-index
# α ∈ ℕᴰ, write the partial sums `A₀ = 0`, `A_k = α₁ + … + α_k`, and
#
#     σ_k = 1 - (x_{k+1} + … + x_D),   so σ₀ = λ₀ = 1 - Σx  and σ_D = 1.
#
# Then, with `Rₙ⁽ᵃ⁾` the *homogeneous* Jacobi polynomial below,
#
#     φ_α(x) = N_α ∏_{k=1..D} R_{α_k}^{(a_k)}(x_k, σ_{k-1}),   a_k = 2A_{k-1} + k-1,
#
# is a polynomial of total degree |α|, and `{φ_α : |α| ≤ K}` is an L²-orthonormal
# basis of `P_K` with `N_α = sqrt( ∏_{k=1..D} (2A_k + k) )`.
#
# This is the Koornwinder construction written on the unit simplex: the usual
# collapsed ("Duffy") coordinates are η_k = x_k/σ_k, and since σ_k = x_k + σ_{k-1}
# the k-th factor is `σ_k^{α_k} P_{α_k}^{(a_k,0)}(2η_k - 1)`. In `D = 1` it
# reduces to `sqrt(2n+1) P_n(2x-1)`, i.e. to `Legendre` exactly.
#
# WHY THE HOMOGENEOUS FORM
#
# Written with η_k the basis has a 0/0 at the vertices where σ_{k-1} = 0 — the
# vertex `e_D` in 2D. The functions themselves are perfectly good polynomials
# there; only the formula is singular. Since the singular point is a vertex,
# quadrature never sees it and the naive form "works" — until an element with
# point-evaluated DoFs asks for a value at a vertex and gets NaN.
#
# So evaluation goes through the homogeneous polynomial
#
#     Rₙ⁽ᵃ⁾(u,v) = (u+v)ⁿ Pₙ⁽ᵃ'⁰⁾((u-v)/(u+v)),
#
# which the Jacobi three-term recurrence turns into, with s = u+v, w = u-v and
# c = 2n+a,
#
#     R₀ = 1,   R₁ = (a+1)u - v,
#     2n(n+a)(c-2) Rₙ = (c-1)[c(c-2)w + a²s] R_{n-1} - 2(n+a-1)(n-1)c s² R_{n-2}.
#
# No division by anything that vanishes on the simplex. Derivatives come from
# differentiating the same recurrence, seeded with `∂ᵤR₁ = a+1`, `∂ᵥR₁ = -1`.
#
# Hessians are not implemented — `PolynomialBasis` allows that, and the intended
# use (moment weights, and prebases used through `linear_combination`) does not
# need them.
############################################################################################

"""
    Dubiner <: Polynomial

The Dubiner (Koornwinder) polynomials, orthonormal for the L² scalar product on
the reference simplex. In 1D they coincide with [`Legendre`](@ref).
"""
struct Dubiner <: Polynomial end

isHierarchical(::Type{Dubiner}) = true

"""
    DubinerBasis{D,V} <: PolynomialBasis{D,V,Dubiner}

Basis of Dubiner polynomials on the reference `D`-simplex, scalar valued or,
for a `MultiValue` `V`, the direct sum of one copy per independent component.

The basis polynomial of multi-index `α` has total degree `|α|`, so the default
`_p_filter` gives `P_order` and any filter selects a subspace by degree.

# Fields
- `order`  the maximum total degree
- `terms`  the multi-indices, as 1-based `CartesianIndex`, as in `CartProdPolyBasis`
- `scales` the L² normalisation `N_α` of each term
"""
struct DubinerBasis{D,V} <: PolynomialBasis{D,V,Dubiner}
  order::Int
  terms::Vector{CartesianIndex{D}}
  scales::Vector{Float64}

  function DubinerBasis{D}(
    ::Type{V}, order::Int, terms::Vector{CartesianIndex{D}}) where {D,V}

    VV = make_concretetype(V)
    scales = [_dubiner_scale(Tuple(t) .- 1) for t in terms]
    new{D,VV}(order, terms, scales)
  end
end

# N_α = sqrt(∏_k (2A_k + k)), the constant making φ_α a unit vector of L² of the
# reference simplex, where A_k = α₁ + … + α_k. It follows from the collapsed-
# coordinate Jacobian ∏_k σ_{k-1} together with
# ∫₀¹ (1-y)ᵃ Pₙ⁽ᵃ'⁰⁾(2y-1)² dy = 1/(2n+a+1); the k-th factor contributes
# 1/(2A_k + k).
function _dubiner_scale(α::NTuple{D,Int}) where D
  s, A = 1.0, 0
  for k in 1:D
    A += α[k]
    s *= (2*A + k)
  end
  return sqrt(s)
end

"""
    DubinerBasis(::Val{D}, ::Type{V}, order::Int [, filter::Function])
    DubinerBasis(::Val{D}, ::Type{V}, order::Int, terms::Vector{CartesianIndex{D}})

Constructors for [`DubinerBasis`](@ref). The filter defaults to `_p_filter`,
i.e. to a basis of `P_order`; it is applied exactly as for any other basis of
the module, through `_define_terms`.
"""
function DubinerBasis(
  ::Val{D}, ::Type{V}, order::Int, filter::Function=_p_filter) where {D,V}

  terms = _define_terms(filter, tfill(order, Val(D)))
  DubinerBasis{D}(V, order, terms)
end

function DubinerBasis(
  ::Val{D}, ::Type{V}, order::Int, terms::Vector{CartesianIndex{D}}) where {D,V}

  DubinerBasis{D}(V, order, terms)
end

Base.size(b::DubinerBasis{D,V}) where {D,V} =
  (length(b.terms) * num_indep_components(V),)

get_order(b::DubinerBasis) = b.order
get_orders(b::DubinerBasis{D}) where D = tfill(get_order(b), Val(D))
get_exponents(b::DubinerBasis) = [Tuple(t) .- 1 for t in b.terms]

testvalue(::Type{DubinerBasis{D,V}}) where {D,V} = DubinerBasis(Val(D), V, 0)

# Caches: one table per spatial direction of the homogeneous Jacobi polynomials
# Rₙ⁽ᵃ⁾(x_d, σ_{d-1}), indexed [d, a+1, n+1]. The Jacobi parameter is
# a_k = 2A_{k-1} + k-1 ≤ 2K + D - 1, so 2K + D values of a suffice.
#
# This is why the default `_return_cache` is overridden: it allocates the
# D × (K+1) table of a tensor-product basis, whose 1D polynomials do not depend
# on the other directions' indices. Here they do, through a.

@inline _dubiner_na(K, D) = 2*K + D

function _return_cache(
  b::DubinerBasis{D}, x, ::Type{G}, ::Val{N_deriv}) where {D,G,N_deriv}

  @notimplementedif N_deriv > 1 """\n
  `DubinerBasis` implements values and gradients, not $(N_deriv)th derivatives.
  """
  T = eltype(G)
  K = get_order(b)
  na = _dubiner_na(K, D)

  r = CachedArray(zeros(G, (length(x), length(b))))
  s = MArray{Tuple{Vararg{D,N_deriv}},T}(undef)
  c = CachedArray(zeros(T, (D, na, K + 1)))
  # the derivative table holds ∂ᵤR and ∂ᵥR in its last index
  t = ntuple(_ -> CachedArray(zeros(T, (D, na, K + 1, 2))), Val(N_deriv))
  return (r, s, c, t...)
end

function _setsize!(b::DubinerBasis{D}, np, r, t...) where D
  K = get_order(b)
  na = _dubiner_na(K, D)
  setsize!(r, (np, length(b)))
  setsize!(t[1], (D, na, K + 1))
  for i in 2:length(t)
    setsize!(t[i], (D, na, K + 1, 2))
  end
end

# The pair (x_d, σ_{d-1}) the d-th factor is a homogeneous Jacobi polynomial of,
# with σ_{d-1} = 1 - (x_d + … + x_D). Note σ₀ is the barycentric coordinate λ₀
# and σ_{D-1} = 1 - x_D, so the last factor's u + v is 1 and it is a plain
# shifted Jacobi polynomial.
@inline function _dubiner_args(x, d, ::Val{D}) where D
  u = x[d]
  v = one(u)
  for i in d:D
    v -= x[i]
  end
  return (u, v)
end

# Fill c[d,a+1,n+1] = Rₙ⁽ᵃ⁾(x_d, σ_{d-1}) for every direction d, every Jacobi
# parameter a in 0:2K+D-1 and every degree n in 0:K, by the three-term
# recurrence given in the header.
function _dubiner_tables!(c::AbstractArray{T,3}, x, ::Val{D}, K) where {T,D}
  na = size(c, 2)
  for d in 1:D
    u, v = _dubiner_args(x, d, Val(D))
    s, w = u + v, u - v
    for a in 0:(na - 1)
      @inbounds c[d, a + 1, 1] = one(T)
      K == 0 && continue
      @inbounds c[d, a + 1, 2] = (a + 1) * u - v
      for n in 2:K
        cc = 2*n + a
        E = 2*n * (n + a) * (cc - 2)
        F = (cc - 1) * (cc * (cc - 2) * w + a * a * s)
        H = 2*(n + a - 1) * (n - 1) * cc
        @inbounds c[d, a + 1, n + 1] =
          (F * c[d, a + 1, n] - H * s * s * c[d, a + 1, n - 1]) / E
      end
    end
  end
end

# As above, and additionally g[d,a+1,n+1,1] = ∂ᵤRₙ⁽ᵃ⁾ and
# g[d,a+1,n+1,2] = ∂ᵥRₙ⁽ᵃ⁾, by differentiating the same recurrence in u and in v
# (∂ᵤw = 1, ∂ᵥw = -1, ∂ᵤs = ∂ᵥs = 1).
function _dubiner_tables!(
  c::AbstractArray{T,3}, g::AbstractArray{T,4}, x, ::Val{D}, K) where {T,D}

  na = size(c, 2)
  z = zero(T)
  for d in 1:D
    u, v = _dubiner_args(x, d, Val(D))
    s, w = u + v, u - v
    for a in 0:(na - 1)
      @inbounds c[d, a + 1, 1] = one(T)
      @inbounds g[d, a + 1, 1, 1] = z
      @inbounds g[d, a + 1, 1, 2] = z
      K == 0 && continue
      @inbounds c[d, a + 1, 2] = (a + 1) * u - v
      @inbounds g[d, a + 1, 2, 1] = a + 1
      @inbounds g[d, a + 1, 2, 2] = -one(T)
      for n in 2:K
        cc = 2*n + a
        E  = 2*n * (n + a) * (cc - 2)
        F  = (cc - 1) * (cc * (cc - 2) * w + a * a * s)
        Fu = (cc - 1) * (cc * (cc - 2) + a * a)
        Fv = (cc - 1) * (a * a - cc * (cc - 2))
        H  = 2*(n + a - 1) * (n - 1) * cc
        @inbounds begin
          cm1, cm2 = c[d, a + 1, n], c[d, a + 1, n - 1]
          c[d, a + 1, n + 1] = (F * cm1 - H * s * s * cm2) / E
          for (q, Fq) in ((1, Fu), (2, Fv))
            g[d, a + 1, n + 1, q] =
              (Fq * cm1 + F * g[d, a + 1, n, q]
               - H * (2*s * cm2 + s * s * g[d, a + 1, n - 1, q])) / E
          end
        end
      end
    end
  end
end

# Gather into F the D factors R_{α_k}^{(a_k)} of the term α (1-based) out of the
# table c, walking the partial sums A_k to get each a_k = 2A_{k-1} + k-1, which
# is left in `as` for the gradient's use.
@inline function _dubiner_factors!(F, as, c, α::NTuple{D,Int}, ::Val{D}) where D
  A = 0
  @inbounds for d in 1:D
    a = 2*A + (d - 1)
    as[d] = a
    F[d] = c[d, a + 1, α[d]]
    A += α[d] - 1
  end
  nothing
end

function _evaluate_nd!(
  b::DubinerBasis{D,V}, x, r::AbstractMatrix, i,
  c::AbstractArray{T,3}, K) where {D,V,T}

  _dubiner_tables!(c, x, Val(D), K)

  F  = MVector{D,T}(undef)
  as = MVector{D,Int}(undef)
  k = 1
  for (it, ci) in enumerate(b.terms)
    _dubiner_factors!(F, as, c, Tuple(ci), Val(D))
    s = T(b.scales[it])
    for d in 1:D
      s *= F[d]
    end
    k = _cartprod_set_value!(r, i, s, k)
  end
end

function _gradient_nd!(
  b::DubinerBasis{D,V}, x, r::AbstractMatrix{G}, i,
  c::AbstractArray{T,3}, g::AbstractArray{T,4},
  s::MVector{D,T}, K) where {D,V,G,T}

  _dubiner_tables!(c, g, x, Val(D), K)

  F  = MVector{D,T}(undef)   # the D factors
  Ru = MVector{D,T}(undef)   # ∂ᵤ of each factor
  Rv = MVector{D,T}(undef)   # ∂ᵥ of each factor
  O  = MVector{D,T}(undef)   # the product of all factors but one
  as = MVector{D,Int}(undef) # the Jacobi parameters

  k = 1
  for (it, ci) in enumerate(b.terms)
    α = Tuple(ci)
    _dubiner_factors!(F, as, c, α, Val(D))
    @inbounds for d in 1:D
      Ru[d] = g[d, as[d] + 1, α[d], 1]
      Rv[d] = g[d, as[d] + 1, α[d], 2]
    end

    # prefix/suffix products, so that a vanishing factor is harmless
    @inbounds begin
      O[1] = one(T)
      for d in 2:D
        O[d] = O[d - 1] * F[d - 1]
      end
      suf = one(T)
      for d in D:-1:1
        O[d] *= suf
        suf *= F[d]
      end
    end

    # ∂_j φ = N Σ_k (∂_j F_k) ∏_{m≠k} F_m, and since u_k = x_k while
    # σ_{k-1} = 1 - (x_k + … + x_D), the k-th factor's derivative is
    #   ∂_j F_k = δ_{jk} ∂ᵤR_k - [j ≥ k] ∂ᵥR_k .
    N = T(b.scales[it])
    @inbounds for j in 1:D
      acc = zero(T)
      for d in 1:D
        dF = ifelse(j == d, Ru[d], zero(T))
        dF = ifelse(j >= d, dF - Rv[d], dF)
        acc += O[d] * dF
      end
      s[j] = N * acc
    end

    k = _cartprod_set_derivative!(r, i, s, k, V)
  end
end

function _hessian_nd!(b::DubinerBasis, x, r, i, c, g, h, s, K)
  @notimplemented "`DubinerBasis` does not implement hessians."
end
