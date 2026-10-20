"""
    SerendipityPolyBasis{D,V,PT} <: PolynomialBasis{D,V,PT}

Basis of the vector proxy (`V<:VectorValue{L}`, `L` = binomial(`D`,`k`)) of the
serendipity space of polynomial `k`-forms on the `D`-cube [AA, (17)]

  𝓢ᵣΛᵏ  = 𝓟ᵣΛᵏ ⊕ 𝓙ᵣΛᵏ ⊕ d𝓙ᵣ₊₁Λᵏ⁻¹,

or of its trimmed counterpart [GK, (3.4)]

  𝓢⁻ᵣΛᵏ = 𝓟ᵣ₋₁Λᵏ ⊕ κ𝓗ᵣ₋₁Λᵏ⁺¹ ⊕ 𝓙ᵣΛᵏ ⊕ d𝓙ᵣΛᵏ⁻¹,

where κ is the Koszul operator, d the exterior derivative, 𝓗ₛ the homogeneous
polynomials of degree s, and 𝓙ᵣΛᵏ is spanned by the κm for the monomial
(`k`+1)-forms m = xᵅdx_σ such that

  deg m ≥ r   and   deg m - ldeg m ≤ r-1,

ldeg m counting the variables xᵢ, i ∉ σ, that enter m linearly [AA, Prop. 3.1].

The basis polynomials are, in this order,

  φα dx_I,    κ(φα dx_σ),    dκ(φα dx_σ),

with φα(x) = φα₁(x₁) × ... × φα`D`(x`D`) a product of 1D polynomials of the
family `PT`. The first kind spans 𝓟ᵣΛᵏ (or 𝓟ᵣ₋₁Λᵏ) like a
[`CartProdPolyBasis`](@ref), the other two span the Koszul parts.

`PT` must be hierarchical: φαᵢ then differs from x^αᵢ by lower powers only, and
these keep each generator in the space, so the index sets defining the space
by monomials define it for any hierarchical family.

The Koszul generators are dependent (κκ = 0). They are made independent by
keeping κ(xᵅdx_σ) only if σ contains the index a₀ = min{i : γᵢ ≥ 2} of the
multidegree γ = α + 1_σ (or a₀ = min{i : γᵢ ≥ 1} if all γᵢ ≤ 1): κ preserves γ
and acts within it as the simplicial boundary on the subsets of supp(γ), whose
image is spanned by the cone over a₀, and choosing γ(a₀) ≥ 2 keeps the cone
inside the subspace selected by the 𝓙 condition. For the whole of 𝓗ₛ the
rule reduces to the one of [AFW] (αᵢ = 0 for i < min σ). The exterior
derivative is injective on the image of κ, so the same rule serves d𝓙.

κ = ι_{x-c} is taken about the centre c = (½,...,½) of the reference cube
[0,1]ᴰ. The spaces are translation invariant, so this changes the basis and not
the space, and it keeps the parity of the 1D families about ½ (better
conditioned Vandermonde matrices).

The forms are transformed into vectors like in [`BarycentricPΛBasis`](@ref),
see `_basis_forms_components`. First order derivatives are implemented.

References:
- [AA] D.N. Arnold & G. Awanou, Finite element differential forms on cubical
  meshes, Math. Comp. 83 (2014)
- [GK] A. Gillette & T. Kloefkorn, Trimmed serendipity finite element
  differential forms, Math. Comp. 88 (2019)
- [AFW] D.N. Arnold, R.S. Falk & R. Winther, Geometric decompositions and local
  bases for spaces of finite element differential forms, CMAME 198 (2009)
"""
struct SerendipityPolyBasis{D,V,PT} <: PolynomialBasis{D,V,PT}
  max_order::Int
  # degree of the 1D polynomials φαᵢ, that is, of the cached 1D tables
  max_order_1d::Int
  # α of the φα dx_I, in index base 1 like the terms of a CartProdPolyBasis
  pterms::Vector{CartesianIndex{D}}
  # (α, σ_id) of the κ(φα dx_σ) and dκ(φα dx_σ), σ_id = combination_index(σ,D)
  κterms::Vector{Tuple{CartesianIndex{D},Int}}
  dκterms::Vector{Tuple{CartesianIndex{D},Int}}
  # κ(dx_σ)    = Σ a (xₗ-½) dx_I,       over (I, a, l)    ∈ κdx[σ_id]
  κdx::Vector{Vector{NTuple{3,Int}}}
  # dκ(f dx_σ) = Σ a (xₗ-½) ∂_q f dx_I, over (I, a, l, q) ∈ dκdx[σ_id],
  #              an entry with l = q = 0 standing for a f dx_I
  dκdx::Vector{Vector{NTuple{4,Int}}}

  function SerendipityPolyBasis{D}(
    ::Type{PT}, ::Type{V}, max_order::Int, max_order_1d::Int,
    pterms::Vector{CartesianIndex{D}},
    κterms::Vector{Tuple{CartesianIndex{D},Int}},
    dκterms::Vector{Tuple{CartesianIndex{D},Int}},
    κdx::Vector{Vector{NTuple{3,Int}}},
    dκdx::Vector{Vector{NTuple{4,Int}}}) where {D,PT<:Polynomial,V}

    @check isconcretetype(PT) "PT needs to be a concrete <:Polynomial type"
    @check isHierarchical(PT) "The polynomial basis must be hierarchical for this space."
    new{D,V,PT}(max_order,max_order_1d,pterms,κterms,dκterms,κdx,dκdx)
  end
end

"""
    SerendipityPolyBasis{D}(PT, T, r, k; trimmed=false, rotate_90=false)

Basis of 𝓢ᵣΛᵏ(□ᴰ), or of 𝓢⁻ᵣΛᵏ(□ᴰ) if `trimmed`, for 1 ≤ `k` ≤ `D`-1 and
`r` ≥ 1, with `T` the scalar component type. `rotate_90` is described in
[`FEEC_poly_basis`](@ref).
"""
function SerendipityPolyBasis{D}(
  ::Type{PT}, ::Type{T}, r::Int, k::Int;
  trimmed=false, rotate_90=false) where {D,PT<:Polynomial,T}

  FEEC_space_definition_checks(Val(D),T,r,k,ifelse(trimmed, :S⁻, :S),rotate_90)
  @check 0 < k < D "Only 1 ≤ k ≤ D-1 is implemented, got k=$k and D=$D"
  @check r ≥ 1 "The polynomial order r must be at least 1, got $r"

  V = VectorValue{binomial(D,k),T}
  proxy = _form_proxy(D,k,rotate_90)

  # 𝓟ᵣΛᵏ, or 𝓟ᵣ₋₁Λᵏ
  pterms = _define_terms(_p_filter, tfill(ifelse(trimmed, r-1, r), Val(D)))

  # 𝓙ᵣΛᵏ, and κ𝓗ᵣ₋₁Λᵏ⁺¹ if trimmed
  κmonomials = _serendipity_J_terms(D,r,k)
  if trimmed
    append!(κmonomials, [ (α, σ) for α in _monomial_exponents(D, r-1)
                                 for σ in sorted_combinations(D,k+1)
                                 if _koszul_cone_filter(α,σ) ])
  end
  κterms = [ (CartesianIndex(α .+ 1), combination_index(σ,D)) for (α, σ) in κmonomials ]
  κdx = [ _koszul_coefficients(σ,proxy) for σ in sorted_combinations(D,k+1) ]

  # d𝓙ᵣ₊₁Λᵏ⁻¹, or d𝓙ᵣΛᵏ⁻¹ if trimmed
  dκmonomials = _serendipity_J_terms(D,ifelse(trimmed, r, r+1),k-1)
  dκterms = [ (CartesianIndex(α .+ 1), combination_index(σ,D)) for (α, σ) in dκmonomials ]
  dκdx = [ _exterior_derivative_koszul_coefficients(D,σ,proxy) for σ in sorted_combinations(D,k) ]

  deg(ci) = sum(Tuple(ci)) - D
  max_order = max(maximum(deg, pterms),
                  maximum(t -> deg(first(t)) + 1, κterms; init=0),
                  maximum(t -> deg(first(t)), dκterms; init=0))
  all_terms = Iterators.flatten((pterms, first.(κterms), first.(dκterms)))
  max_order_1d = maximum(ci -> maximum(Tuple(ci)) - 1, all_terms)

  SerendipityPolyBasis{D}(PT,V,max_order,max_order_1d,pterms,κterms,dκterms,κdx,dκdx)
end

function Base.size(b::SerendipityPolyBasis{D,V}) where {D,V}
  (num_indep_components(V)*length(b.pterms) + length(b.κterms) + length(b.dκterms),)
end
get_order(b::SerendipityPolyBasis) = b.max_order
get_orders(b::SerendipityPolyBasis{D}) where D = tfill(get_order(b), Val(D))
_get_parameters(b::SerendipityPolyBasis) = b.max_order_1d

function testvalue(::Type{SerendipityPolyBasis{D,V,PT}}) where {D,V,PT}
  terms = Tuple{CartesianIndex{D},Int}[]
  SerendipityPolyBasis{D}(PT,V,0,0,CartesianIndex{D}[],terms,copy(terms),
                          Vector{NTuple{3,Int}}[],Vector{NTuple{4,Int}}[])
end

"""
    _exact_poly_basis(PT, T, D, s, k; rotate_90=false)

The generators dκ(φα dx_σ), |α| = `s`, of a [`SerendipityPolyBasis`](@ref): a
basis of the vector proxy of d𝓗ₛ₊₁Λᵏ⁻¹(ℝᴰ) = dκ𝓗ₛΛᵏ, up to terms in 𝓟ₛ₋₁Λᵏ.
"""
function _exact_poly_basis(::Type{PT}, ::Type{T}, D, s, k; rotate_90=false) where {PT,T}
  V = VectorValue{binomial(D,k),T}
  proxy = _form_proxy(D,k,rotate_90)
  dκterms = [ (CartesianIndex(α .+ 1), combination_index(σ,D))
              for α in _monomial_exponents(D, s) for σ in sorted_combinations(D,k)
              if _koszul_cone_filter(α,σ) ]
  dκdx = [ _exterior_derivative_koszul_coefficients(D,σ,proxy) for σ in sorted_combinations(D,k) ]
  SerendipityPolyBasis{D}(PT,V,s,s,CartesianIndex{D}[],Tuple{CartesianIndex{D},Int}[],
                          dκterms,Vector{NTuple{3,Int}}[],dκdx)
end

###########################
# Index sets of the space #
###########################

# The multi-indices α ∈ ℕᴰ with |α| ∈ degrees
function _monomial_exponents(D, degrees)
  s = maximum(degrees)
  exponents = (Tuple(ci) .- 1 for ci in CartesianIndices(tfill(s+1, Val(D))))
  [ α for α in exponents if sum(α) in degrees ]
end

# Number of variables xᵢ, i ∉ σ, entering xᵅdx_σ linearly
_ldeg(α, σ) = count(i -> !(i in σ) && isone(α[i]), eachindex(α))

"""
    _koszul_cone_filter(α, σ)

Whether κ(xᵅdx_σ) is kept in a basis of the Koszul image, that is, whether
σ ∋ a₀ = min{i : γᵢ ≥ 2}, or min{i : γᵢ ≥ 1} if there is none, for the
multidegree γ = α + 1_σ, see [`SerendipityPolyBasis`](@ref).
"""
function _koszul_cone_filter(α, σ)
  γ = ntuple(i -> α[i] + Int(i in σ), length(α))
  a₀ = something(findfirst(≥(2), γ), findfirst(≥(1), γ))
  a₀ in σ
end

"""
    _serendipity_J_terms(D, r, k)

The monomial (`k`+1)-forms xᵅdx_σ, as (α, σ), whose Koszul images form a basis
of 𝓙ᵣΛᵏ(ℝᴰ), see [`SerendipityPolyBasis`](@ref).
"""
function _serendipity_J_terms(D, r, k)
  # ldeg ≤ D-k-1, so the 𝓙 condition bounds the degree by r+D-k-1
  [ (α, σ) for α in _monomial_exponents(D, r:(r+D-k-1))
           for σ in sorted_combinations(D,k+1)
           if sum(α) - _ldeg(α,σ) ≤ r-1 && _koszul_cone_filter(α,σ) ]
end

##########################################
# Coefficients of the Koszul basis forms #
##########################################

# I => (proxy component, sign) for the k-combinations I, see _basis_forms_components
function _form_proxy(D, k, rotate_90)
  components = _basis_forms_components(D,k,false,rotate_90)
  Dict(I => (I_id, sgn) for (I_id, I, sgn) in components)
end

# κ(dx_σ) = Σⱼ (-1)ʲ⁻¹ (x_σⱼ-½) dx_{σ∖σⱼ}
function _koszul_coefficients(σ, proxy)
  map(eachindex(σ)) do j
    I_id, sgn = proxy[deleteat!(copy(σ),j)]
    (I_id, (-1)^(j-1)*sgn, σ[j])
  end
end

# dκ(f dx_σ) = k f dx_σ + df ∧ κ(dx_σ), since dκ(dx_σ) = k dx_σ, so
#            = k f dx_σ + Σⱼ Σ_q (-1)ʲ⁻¹ (x_σⱼ-½) ∂_q f dx_q ∧ dx_{σ∖σⱼ}
function _exterior_derivative_koszul_coefficients(D, σ, proxy)
  σ_id, σ_sgn = proxy[σ]
  coeffs = [(σ_id, length(σ)*σ_sgn, 0, 0)]
  for j in eachindex(σ)
    τ = deleteat!(copy(σ),j)
    for q in 1:D
      q in τ && continue
      I_id, sgn = proxy[sort!([q; τ])]
      push!(coeffs, (I_id, (-1)^(j-1)*sorting_sign(q, τ...)*sgn, σ[j], q))
    end
  end
  coeffs
end

#################################
# nD evaluations implementation #
#################################

# The dκ generators differentiate the 1D polynomials once more than the basis
# itself is differentiated, so N_deriv+2 tables of 1D derivatives are cached
# and passed together as `c`. For gradients, `g` is a mutable accumulator for
# the gradient of one basis polynomial.
function _return_cache(
  b::SerendipityPolyBasis{D}, x, ::Type{G}, ::Val{N_deriv}) where {D,G,N_deriv}

  # third 1D derivatives would be needed, which no Polynomial provides
  @notimplementedif N_deriv > 1 "Hessians of SerendipityPolyBasis are not implemented"
  T = eltype(G)
  np = length(x)
  ndof = length(b)
  ndof_1d = b.max_order_1d + 1
  r = CachedArray(zeros(G,(np,ndof)))
  s = MArray{Tuple{Vararg{D,N_deriv}},T}(undef)
  t = ntuple( _ -> CachedArray(zeros(T,(D,ndof_1d))), Val(N_deriv+2))
  iszero(N_deriv) ? (r, s, t) : (r, s, t, zero(Mutable(G)))
end

function _setsize!(b::SerendipityPolyBasis{D}, np, r, t, _...) where D
  setsize!(r,(np,length(b)))
  for c in t
    setsize!(c,(D,b.max_order_1d+1))
  end
end

# The values at x of the D 1D factors of φα, read from a table c of 1D
# polynomials (or of their derivatives).
@inline function _factors(c, ci::CartesianIndex{D}) where D
  ntuple(d -> @inbounds(c[d,ci[d]]), Val(D))
end

# ∂_q φα and ∂_p ∂_q φα, from the values v, derivatives g and second
# derivatives h of the 1D factors of φα.
@inline function _∂(v::NTuple{D}, g, q) where D
  prod(ntuple(d -> ifelse(d == q, g[d], v[d]), Val(D)))
end
@inline function _∂∂(v::NTuple{D}, g, h, p, q) where D
  prod(ntuple(d -> ifelse(d == p == q, h[d], ifelse(d == p || d == q, g[d], v[d])), Val(D)))
end

function _evaluate_nd!(
  b::SerendipityPolyBasis{D,V,PT}, x,
  r::AbstractMatrix{Vr}, i,
  t::NTuple{2}, K) where {D,V,PT,Vr}

  for d in 1:D
    _derivatives_1d!(PT,K,t,x,d)
  end
  c, g = t
  T = eltype(Vr)
  ξ = Tuple(x) .- T(1/2) # κ is taken about the cube centre
  m = zero(Mutable(Vr))
  k = 1

  @inbounds begin
    for ci in b.pterms
      φ = prod(_factors(c,ci))
      k = _cartprod_set_value!(r,i,φ,k)
    end

    for (ci, σ) in b.κterms
      φ = prod(_factors(c,ci))
      m .= zero(T)
      for (I, a, l) in b.κdx[σ]
        m[I] += a * ξ[l] * φ
      end
      r[i,k] = m
      k += 1
    end

    for (ci, σ) in b.dκterms
      v, ∂v = _factors(c,ci), _factors(g,ci)
      m .= zero(T)
      for (I, a, l, q) in b.dκdx[σ]
        m[I] += a * (iszero(l) ? prod(v) : ξ[l] * _∂(v,∂v,q))
      end
      r[i,k] = m
      k += 1
    end
  end
end

function _gradient_nd!(
  b::SerendipityPolyBasis{D,V,PT}, x,
  r::AbstractMatrix{G}, i,
  t::NTuple{3}, m::MArray,
  s::MVector{D,T}, K) where {D,V,PT,G,T}

  for d in 1:D
    _derivatives_1d!(PT,K,t,x,d)
  end
  c, g, h = t
  ξ = Tuple(x) .- T(1/2)
  k = 1

  @inbounds begin
    for ci in b.pterms
      v, ∂v = _factors(c,ci), _factors(g,ci)
      for q in 1:D
        s[q] = _∂(v,∂v,q)
      end
      k = _cartprod_set_derivative!(r,i,s,k,V)
    end

    # ∂_p((xₗ-½) φα) = δₚₗ φα + (xₗ-½) ∂_p φα
    for (ci, σ) in b.κterms
      v, ∂v = _factors(c,ci), _factors(g,ci)
      φ = prod(v)
      m .= zero(T)
      for (I, a, l) in b.κdx[σ], p in 1:D
        m[p,I] += a * (ξ[l] * _∂(v,∂v,p) + ifelse(p == l, φ, zero(T)))
      end
      r[i,k] = m
      k += 1
    end

    # ∂_p((xₗ-½) ∂_q φα) = δₚₗ ∂_q φα + (xₗ-½) ∂_p ∂_q φα
    for (ci, σ) in b.dκterms
      v, ∂v, ∂²v = _factors(c,ci), _factors(g,ci), _factors(h,ci)
      m .= zero(T)
      for (I, a, l, q) in b.dκdx[σ], p in 1:D
        if iszero(l)
          m[p,I] += a * _∂(v,∂v,p)
        else
          m[p,I] += a * (ξ[l] * _∂∂(v,∂v,∂²v,p,q) + ifelse(p == l, _∂(v,∂v,q), zero(T)))
        end
      end
      r[i,k] = m
      k += 1
    end
  end
end
