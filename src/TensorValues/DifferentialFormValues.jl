# DifferentialFormValues.jl
#
# Defines DifferentialFormValue together with its pointwise algebra
# (∧, ι, ⋆, ♭/♯, koszul, …).
#
# Operations needing a derivative are not part of this pointwise algebra.
# `exterior_derivative` and `codifferential` belong to Gridap.Fields; the
# generic functions declared here (`lie_derivative`, `symbolic_coordinates`)
# get their methods from the GridapSymbolicsExt package extension, which loads
# when Symbolics is present.

using Combinatorics: levicivita

"""
    DifferentialFormValue{K,D,T,L} <: MultiValue{NTuple{K,D},T,K,L}

Value of a differential K-form in D dimensions: the `L = binomial(D,K)`
components on the orientation-ordered basis `{dx^I}` (lexicographic
K-combinations `I` of `1:D`), with scalar type `T`.

The components are coefficients on a coframe; which coframe they refer to is a
property of the space the value came from, not of the value. `show` labels them
`dxⁱ`, or `dλⁱ` when the `IOContext` property `:coordinates` is `:barycentric`.
"""
struct DifferentialFormValue{K,D,T,L} <: MultiValue{NTuple{K,D},T,K,L}
  data::NTuple{L,T}

  # All-params inner constructor: no validation.
  function DifferentialFormValue{K,D,T,L}(data::NTuple{L,T}) where {K,D,T,L}
    new{K,D,T,L}(data)
  end

  # 2-param shorthand (validates component count).
  function DifferentialFormValue{K,D}(data::NTuple{L,T}) where {K,D,T,L}
    @assert L == binomial(D,K) "wrong number of values: got $L, expected $(binomial(D,K))"
    new{K,D,T,L}(data)
  end

  # Empty tuple (K > D): no component.
  function DifferentialFormValue{K,D}(data::Tuple{}) where {K,D}
    @assert binomial(D,K) == 0 "Empty tuple requires K > D"
    new{K,D,Float64,0}(data)
  end
end

# ============================================================
# Display
# ============================================================

const _sbs_cart = ["dx¹","dx²","dx³","dx⁴","dx⁵","dx⁶","dx⁷","dx⁸","dx⁹"]
const _sbs_bary = ["dλ¹","dλ²","dλ³","dλ⁴","dλ⁵","dλ⁶","dλ⁷","dλ⁸","dλ⁹"]

function _show_dfv(io::IO, a::DifferentialFormValue{K,D}, basis_strs) where {K,D}
  L = length(a.data)
  L == 0 && return print(io, "0")
  cs = sorted_combinations(D, K)
  parts = [join(basis_strs[I], " ∧ ") for I in cs]
  result = [string(ai, " ", bi) for (ai, bi) in zip(a.data, parts)]
  print(io, join(result, " + "))
end

function _coframe_labels(io::IO)
  c = get(io, :coordinates, :cartesian)
  c === :cartesian   && return _sbs_cart
  c === :barycentric && return _sbs_bary
  throw(ArgumentError(
    "unknown :coordinates value $(repr(c)); expected :cartesian or :barycentric"))
end

"""
    show(io::IO, ::MIME"text/plain", ω::DifferentialFormValue)

Print `ω` as a linear combination of coframe elements, labelled `dx¹,…,dx^D`.

Set the `:coordinates` `IOContext` property to `:barycentric` to label them
`dλ¹,…,dλ^D` instead, for a form whose components are coefficients on the
ambient barycentric coframe of a (D−1)-simplex:

    show(IOContext(stdout, :coordinates => :barycentric), MIME("text/plain"), ω)
"""
function Base.show(io::IO, ::MIME"text/plain", a::DifferentialFormValue{K,D}) where {K,D}
  @assert D <= 9 "show not implemented for D > 9"
  _show_dfv(io, a, _coframe_labels(io))
end

# ============================================================
# Basic algebra: zero, +, -, scalar *
# ============================================================

function Base.zero(::Type{DifferentialFormValue{K,D,T,L}}) where {K,D,T,L}
  DifferentialFormValue{K,D,T,L}(ntuple(_ -> zero(T), L))
end

# change_eltype's generic <:Number fallback (Operations.jl) would otherwise
# collapse DifferentialFormValue to its bare scalar type T2, since MultiValue
# <: Number. PolynomialBasis's generic return_cache (_return_val_eltype)
# relies on change_eltype to reconstruct the value type, so this is needed
# for any DifferentialFormValue-valued PolynomialBasis to evaluate correctly.
change_eltype(::Type{<:DifferentialFormValue{K,D,T1,L}}, ::Type{T2}) where {K,D,T1,T2,L} =
  DifferentialFormValue{K,D,T2,L}

Base.zero(ω::DifferentialFormValue) = zero(typeof(ω))

function Base.:+(a::DifferentialFormValue{K,D,T1,L}, b::DifferentialFormValue{K,D,T2,L}) where {K,D,T1,T2,L}
  d = map(+, a.data, b.data)
  DifferentialFormValue{K,D,eltype(d),L}(d)
end

function Base.:-(a::DifferentialFormValue{K,D,T1,L}, b::DifferentialFormValue{K,D,T2,L}) where {K,D,T1,T2,L}
  d = map(-, a.data, b.data)
  DifferentialFormValue{K,D,eltype(d),L}(d)
end

Base.:-(a::DifferentialFormValue{K,D,T,L}) where {K,D,T,L} =
  DifferentialFormValue{K,D,T,L}(map(-, a.data))

function Base.:*(s::Union{Real,Complex}, ω::DifferentialFormValue{K,D,T,L}) where {K,D,T,L}
  d = map(x -> s*x, ω.data)
  DifferentialFormValue{K,D,eltype(d),L}(d)
end

Base.:*(ω::DifferentialFormValue, s::Union{Real,Complex}) = s * ω

# ============================================================
# Helper: multi-index ↔ linear index map
#
# Returns an array `arr` such that arr[i₁,...,iₖ] = position of
# the sorted K-combination (i₁,...,iₖ) in the lexicographic list
# of all K-subsets of {1,...,D}.  For K=0 returns a 0-dim array
# containing 1 (there is exactly one empty combination).
# ============================================================

function _ijk_l(K::Int, D::Int)
  cm  = sorted_combinations(D, K)
  _d  = ntuple(_ -> D, K)
  arr = zeros(Int, _d)
  for (c, idx) in enumerate(cm)
    arr[idx...] = c
  end
  arr
end

# ============================================================
# Exterior (wedge) product ∧
# ============================================================

"""
    ∧(a::DifferentialFormValue{K1,D}, b::DifferentialFormValue{K2,D})

Pointwise exterior (wedge) product, a `DifferentialFormValue{K1+K2,D}`.
"""
function ∧(a::DifferentialFormValue{K1,D,T1,L1}, b::DifferentialFormValue{K2,D,T2,L2}) where {K1,K2,D,T1,T2,L1,L2}
  K  = K1 + K2
  T  = typeof(zero(T1) * zero(T2))
  L  = binomial(D, K)   # = 0 when K > D

  d  = zeros(T, max(L, 1))   # avoid zero-length zeros() for accumulation
  c1 = sorted_combinations(D, K1)
  c2 = sorted_combinations(D, K2)

  if L > 0
    ijk_l1  = _ijk_l(K1, D)
    ijk_l2  = _ijk_l(K2, D)
    ijk_l12 = _ijk_l(K,  D)
    for (n1, i1) in enumerate(c1)
      for (n2, i2) in enumerate(c2)
        i12 = [i1..., i2...]
        l   = ijk_l12[sort(i12)...]
        if l > 0
          d[l] += levicivita(sortperm(i12)) * a.data[ijk_l1[i1...]] * b.data[ijk_l2[i2...]]
        end
      end
    end
  end

  DifferentialFormValue{K,D,T,L}(Tuple(d[1:L]))
end

# ============================================================
# Interior product ι_v(ω): contracts ω(v, -)
# ι_v : Ωᴷ(D) → Ωᴷ⁻¹(D),   (ι_v ω)(w₂,…,wₖ) = ω(v, w₂,…,wₖ)
# ============================================================

"""
    interior_product(v::VectorValue{D}, ω::DifferentialFormValue{K,D})
    ι(v, ω)

Interior product (contraction) `ι_v ω`, a `DifferentialFormValue{K-1,D}`:
`(ι_v ω)(w₂,…,wₖ) = ω(v, w₂,…,wₖ)`.
"""
function interior_product(v::VectorValue{D,Tv}, ω::DifferentialFormValue{K,D,Tw,Lw}) where {K,D,Tv,Tw,Lw}
  @assert K >= 1 "interior product requires K ≥ 1"
  T    = typeof(zero(Tv) * zero(Tw))
  Km1  = K - 1
  L    = binomial(D, Km1)
  d    = zeros(T, max(L, 1))

  c_in      = sorted_combinations(D, K)
  ijk_l_out = _ijk_l(Km1, D)

  for (n, c) in enumerate(c_in)
    for j in 1:K
      idx = c[j]
      rem = [c[i] for i in 1:K if i != j]
      sort!(rem)
      l   = Km1 == 0 ? 1 : ijk_l_out[rem...]
      if l > 0
        sgn = iseven(j - 1) ? 1 : -1
        d[l] += sgn * v[idx] * ω.data[n]
      end
    end
  end

  DifferentialFormValue{Km1,D,T,L}(Tuple(d[1:L]))
end

const ι = interior_product

# ============================================================
# Gradient support: outer(x, ω) for 1-forms.
#
# Gridap derives the gradient value type of a V-valued basis as
# typeof(outer(zero(x), zero(V))) (Fields.gradient_type). Defining outer of a
# point with a 1-form value as the TensorValue{D,D} with entries
# [i,j] = x_i·ω_j makes ∇(basis) follow Gridap's convention ∇u[i,j] = ∂u_j/∂xⁱ
# for the vector proxy of the 1-form.
# ============================================================

function outer(a::VectorValue{D,Ta}, b::DifferentialFormValue{1,D,Tb,L}) where {D,Ta,Tb,L}
  T = promote_type(Ta, Tb)
  TensorValue{D,D,T}(ntuple(k -> a[(k-1)%D+1] * b.data[(k-1)÷D+1], Val(D*D)))
end

# ============================================================
# Hodge star ⋆ (flat Euclidean metric, g = Iᴅ, so √det g = 1)
#
# (⋆ω)_J = ∑_I  ε_{IJ}  ω_I
# where I ranges over K-combinations, J over (D-K)-combinations,
# and ε_{IJ} = levicivita of the permutation (I,J) → sorted (1…D).
# ============================================================

"""
    _hodge_star_matrix(K, D) -> Matrix{Int}

Levi-Civita sign matrix of the flat Hodge star: entry `[m,n] = ε_{I_n J_m}`
for `I_n` the n-th K-combination and `J_m` the m-th (D−K)-combination of
`1:D` (zero when they overlap). Shared by `hodge_star` (here) and
`CodifferentialForm` (Gridap.Fields).
"""
function _hodge_star_matrix(K::Int, D::Int)
  Kc = D - K
  c_in  = sorted_combinations(D, K)
  c_out = sorted_combinations(D, Kc)
  M = zeros(Int, length(c_out), length(c_in))
  for (m, J) in enumerate(c_out), (n, I) in enumerate(c_in)
    perm = [I..., J...]
    length(unique(perm)) == D || continue
    M[m, n] = levicivita(sortperm(perm))
  end
  M
end

"""
    hodge_star(ω::DifferentialFormValue{K,D})
    hodge_star(ω::DifferentialFormValue{K,D}, g_inv::SymTensorValue{D}, sqrt_det_g)
    ⋆(ω)
    ⋆(ω, g_inv, sqrt_det_g)

Hodge star `⋆ω`, a `DifferentialFormValue{D-K,D}`. The one-argument form uses
the flat Euclidean metric; the three-argument form takes a pointwise inverse
metric tensor and `√det(g)`.
"""
function hodge_star(ω::DifferentialFormValue{K,D,T,L}) where {K,D,T,L}
  Kc = D - K
  Lc = binomial(D, Kc)
  M  = _hodge_star_matrix(K, D)
  d  = ntuple(m -> sum(M[m, n] * ω.data[n] for n in 1:L), Lc)
  DifferentialFormValue{Kc,D,T,Lc}(d)
end

const ⋆ = hodge_star

# ============================================================
# Generic functions extended elsewhere
#
# Methods for Symbolics.Num coefficients are provided by the GridapSymbolicsExt
# package extension.
# ============================================================

"""
    lie_derivative(v, ω)
    𝓛(v,ω)

Lie derivative ``𝓛_v ω = d(ι_v ω) + ι_v(dω)`` (Cartan's magic formula).

Requires symbolic (`Symbolics.Num`) coefficients: methods are provided by the
GridapSymbolicsExt package extension when Symbolics is loaded.
"""
function lie_derivative end
const 𝓛 = lie_derivative

"""
    symbolic_coordinates(D::Integer)

The first `D` global Cartesian symbolic coordinate variables `x¹,…,x⁹`
(at most 9), as a tuple of `Symbolics.Num`. These are the variables
`exterior_derivative` differentiates with respect to for symbolic forms.

Provided by the GridapSymbolicsExt package extension when Symbolics is loaded.
"""
function symbolic_coordinates end

# ============================================================
# Riemannian Hodge star  ⋆_g : Ωᴷ → Ωᴰ⁻ᴷ
#
# Algorithm (coordinate formula):
#   Step 1 — raise indices:
#     ω^I = Σ_{J sorted K-index} det([g_inv[I[a],J[b]]]_{a,b}) ω_J
#   Step 2 — Levi-Civita contraction:
#     (⋆_g ω)_L = √det(g) Σ_{I sorted} ε(I,L) ω^I
#
# Passing g_inv (the pointwise inverse metric tensor) and sqrt_det_g
# (the pointwise square root of det g) keeps this purely algebraic.
# ============================================================

# Determinant of the K×K submatrix of a symmetric tensor indexed by I (rows) and J (cols)
function _minor_k(g_inv, I, J, K)
  K == 0 && return one(eltype(g_inv.data))
  K == 1 && return g_inv[I[1], J[1]]
  K == 2 && return g_inv[I[1],J[1]]*g_inv[I[2],J[2]] - g_inv[I[1],J[2]]*g_inv[I[2],J[1]]
  M = [g_inv[I[a], J[b]] for a in 1:K, b in 1:K]
  det(M)
end

function hodge_star(ω::DifferentialFormValue{K,D,T}, g_inv::SymTensorValue{D}, sqrt_det_g) where {K,D,T}
  Kc = D - K
  L  = binomial(D, K)
  Lc = binomial(D, Kc)
  Tout = promote_type(T, typeof(sqrt_det_g), eltype(g_inv.data))

  c_K  = sorted_combinations(D, K)
  c_Kc = sorted_combinations(D, Kc)

  # Step 1: raise K indices of ω
  omega_raised = zeros(Tout, max(L, 1))
  for (n, I) in enumerate(c_K)
    for (m, J) in enumerate(c_K)
      omega_raised[n] += _minor_k(g_inv, I, J, K) * ω.data[m]
    end
  end

  # Step 2: Levi-Civita contraction
  d = zeros(Tout, max(Lc, 1))
  ijk_l_out = _ijk_l(Kc, D)

  for (n, I) in enumerate(c_K)
    for J in c_Kc
      perm = [I..., J...]
      length(unique(perm)) == D || continue
      idx = Kc == 0 ? 1 : ijk_l_out[J...]   # J already sorted
      d[idx] += levicivita(sortperm(perm)) * sqrt_det_g * omega_raised[n]
    end
  end

  DifferentialFormValue{Kc,D}(Tuple(d[1:Lc]))
end

# ============================================================
# Musical isomorphisms  ♭ (flat) and ♯ (sharp)
#
# flat  (♭):  VectorValue{D} → DifferentialFormValue{1,D}
#   (v♭)_i = Σ_j g_{ij} v^j           (lowers an index)
#
# sharp (♯):  DifferentialFormValue{1,D} → VectorValue{D}
#   (ω♯)^i = Σ_j g^{ij} ω_j           (raises an index)
# ============================================================

"""
    flat(v::VectorValue{D}, g::SymTensorValue{D})
    ♭(v, g)

Musical isomorphism ♭: lower the index of the vector `v` with the metric `g`,
giving a `DifferentialFormValue{1,D}`.
"""
function flat(v::VectorValue{D,T}, g::SymTensorValue{D}) where {D,T}
  Tout = promote_type(T, eltype(g.data))
  d = ntuple(i -> sum(g[i,j] * v[j] for j in 1:D), D)
  DifferentialFormValue{1,D}(Tuple{Vararg{Tout,D}}(d))
end

const ♭ = flat

"""
    sharp(ω::DifferentialFormValue{1,D}, g_inv::SymTensorValue{D})
    ♯(ω, g_inv)

Musical isomorphism `♯`: raise the index of the 1-form `ω` with the inverse
metric `g_inv`, giving a `VectorValue{D}`.
"""
function sharp(ω::DifferentialFormValue{1,D,T}, g_inv::SymTensorValue{D}) where {D,T}
  Tout = promote_type(T, eltype(g_inv.data))
  d = ntuple(i -> sum(g_inv[i,j] * ω.data[j] for j in 1:D), D)
  VectorValue{D,Tout}(d)
end

const ♯ = sharp

"""
    apply_form(ω::DifferentialFormValue{K,D,T}, v₁, …, vₖ) → scalar

Evaluate the K-form ω on K tangent vectors v₁,…,vₖ ∈ ℝᴰ:
  ω(v₁,…,vₖ) = Σ_{|I|=K} ω_I · det([vₐ[I[b]]]_{a,b=1}^K)
Each vector must be indexable (Tuple, Vector, VectorValue).
For K=0: no vectors needed; returns `ω.data[1]`.
"""
function apply_form(ω::DifferentialFormValue{K,D,T}, vs...) where {K,D,T}
  @assert length(vs) == K "K=$K form requires K vectors, got $(length(vs))"
  K == 0 && return ω.data[1]
  cs = sorted_combinations(D, K)
  result = zero(promote_type(T, Float64))
  for (idx, I) in enumerate(cs)
    M = [vs[a][I[b]] for a in 1:K, b in 1:K]
    result += ω.data[idx] * det(M)
  end
  result
end

# ============================================================
# Pullback  φ*ω : Ωᴷ(Dn) → Ωᴷ(Dm)
#
# Given the Jacobian J = ∇φ  (TensorValue{Dn,Dm}, J[i,j] = ∂φⁱ/∂ξʲ) of
# a map φ : ℝᴰᵐ → ℝᴰⁿ and a K-form ω on ℝᴰⁿ:
#
#   (φ*ω)_J(ξ) = Σ_I  ω_I(φ(ξ)) · det(J[I, J])
#
# where I is a sorted K-multi-index into {1:Dn} (ambient) and
# J is a sorted K-multi-index into {1:Dm} (chart).
# ============================================================

"""
    pullback(ω::DifferentialFormValue{K,Dn}, J::TensorValue{Dn,Dm})

Pointwise pullback `φ^*ω` of a K-form under a map with Jacobian `J = ∇φ`,
giving a `DifferentialFormValue{K,Dm}`. A Field-level method for lazy
pullbacks is provided in `Gridap.Fields`.
"""
function pullback(ω::DifferentialFormValue{K,Dn,T}, J::TensorValue{Dn,Dm,T2}) where {K,Dn,Dm,T,T2}
  Tout = promote_type(T, T2)
  L    = binomial(Dm, K)
  d    = zeros(Tout, max(L, 1))

  c_I  = sorted_combinations(Dn, K)   # K-combos in ambient space
  c_J  = sorted_combinations(Dm, K)   # K-combos in chart space

  for (m, Jidx) in enumerate(c_J)
    for (n, I) in enumerate(c_I)
      d[m] += _minor_k(J, I, Jidx, K) * ω.data[n]
    end
  end

  DifferentialFormValue{K,Dm}(Tuple(d[1:L]))
end

# ============================================================
# Koszul operator  κ_x(ω) = ι_x(ω)  (value level)
#
# κ : Λᴷ → Λᴷ⁻¹,  κ_x(ω) = interior_product(x, ω)
#
# Homotopy identity on Pᵣ Λᴷ: dκ + κd = (r+K) id
# This is the fundamental algebraic fact behind trimmed polynomial
# form spaces: P⁻ᵣΛᴷ = Pᵣ₋₁Λᴷ ⊕ κ(Hᵣ₋₁Λᴷ⁺¹)
# ============================================================

"""
    koszul(x::VectorValue{D}, ω::DifferentialFormValue{K,D})

Koszul differential ``κ_x(ω) = ι_x(ω)``, the interior product with the position
vector `x`. A Field-level method is provided in `Gridap.Fields`.
"""
koszul(x::VectorValue{D}, ω::DifferentialFormValue{K,D}) where {K,D} =
  interior_product(x, ω)

# ============================================================
# Volume form scalar extraction
#
# vol_coeff : Λᴰ → ℝ
# A top-degree D-form in D dimensions has exactly one component.
# Extracting it gives the scalar integrand for ∫(ω ∧ ⋆η) dΩ.
# ============================================================

"""
    vol_coeff(ω::DifferentialFormValue{D,D})

The single scalar coefficient of a top-degree D-form in D dimensions.
Enables using Gridap's standard measure for integration of differential forms:
    ∫(vol_coeff(ω ∧ ⋆η)) * dΩ
"""
vol_coeff(ω::DifferentialFormValue{D,D,T}) where {D,T} = ω[1]

# ============================================================
# Isomorphism: VectorValue{D} ↔ DifferentialFormValue{1,D}
#
# In a flat Euclidean space, 1-forms and vectors are identified
# via the standard inner product.  These conversions let you
# take a VectorValue FEFunction (e.g., from a Nédélec FESpace)
# and treat it as a genuine differential 1-form for subsequent
# operations (d, ⋆, ∧, …).
# ============================================================

"""
    to_1form(v::VectorValue{D})

Conversion of `v` to a 1-form, a `DifferentialFormValue{1,D}`, assuming flat space
(the metric tensor is the iddentity matrix).
"""
to_1form(v::VectorValue{D,T}) where {D,T} = DifferentialFormValue{1,D}(Tuple(v))

"""
    from_1form(ω::DifferentialFormValue{1,D,T})

Conversion of the 1-form value `ω` to a `VectorValue{D}`, assuming flat space
(the metric tensor is the iddentity matrix).
"""
from_1form(ω::DifferentialFormValue{1,D,T}) where {D,T} = VectorValue{D,T}(Tuple(ω))

# General K-form ↔ VectorValue{binomial(D,K)} (component-wise isomorphism)

"""
    to_Kform(v::VectorValue{L}, ::Val{K}, ::Val{D})

Reinterpret the `L = binomial(D,K)` components of `v` as the components of a
`DifferentialFormValue{K,D}` in canonical Cartesian basis `dx^I`.
"""
to_Kform(v::VectorValue{L,T}, ::Val{K}, ::Val{D}) where {L,T,K,D} = DifferentialFormValue{K,D}(Tuple(v))

"""
    from_Kform(ω::DifferentialFormValue{K,D})

Reinterpret the components of `ω` as a `VectorValue{binomial(D,K)}`, the inverse
of [`to_Kform`](@ref).
"""
from_Kform(ω::DifferentialFormValue{K,D,T}) where {K,D,T} = VectorValue{binomial(D,K),T}(Tuple(ω))

# ============================================================
# Scalar ↔ 0-form and D-form
# ============================================================

"""
    to_0form(u::Number, ::Val{D})

Convert the scalar `u` to a `DifferentialFormValue{0,D}`.
"""
to_0form(u::_Scalar, ::Val{D}) where D = DifferentialFormValue{0,D}((u,))

"""
    to_Dform(u::Number, ::Val{D})

Convert the scalar `u` as the sole component of the top-degree
`DifferentialFormValue{D,D}`, the inverse of [`vol_coeff`](@ref).
"""
to_Dform(u::_Scalar, ::Val{D}) where D = DifferentialFormValue{D,D}((u,))

# ============================================================
# L2 inner product of K-forms (flat Euclidean metric)
#
# ⟨ω, η⟩ = Σ_{|I|=K} ω_I η_I   (sum of component products)
#
# This is the isomorphism Λᴷ ≅ ℝ^{binomial(D,K)} as inner product spaces,
# using the standard orientation-ordered basis {dx^I}.
#
# MomentBasedDofBasis uses `⋅` (= LinearAlgebra.dot) to contract moments
# with prebasis values: moment_I ⋅ value_J → Float64 (DOF matrix entry).
# ============================================================

LinearAlgebra.dot(a::DifferentialFormValue{K,D,T,L},
                  b::DifferentialFormValue{K,D,T,L}) where {K,D,T,L} =
  sum(a.data[i] * b.data[i] for i in 1:L)

# Gridap uses ⊙ (inner, full tensor contraction) in MomentBasedDofBasis for DOF evaluation:
#   T = typeof(zero(V) ⊙ zero(Vr))   and   dofs[o] += moments[i,j] ⊙ vals[nodes[i]]
#
# For K-forms, ⊙ must give the flat Euclidean inner product: Σ_I ω_I η_I → Float64.
# Without this override, the generic contracted_product(Val{K}, ...) is called, which
# fails for K≥2 because DifferentialFormValue stores only binomial(D,K)
# antisymmetric components, not the full D^K tensor that contracted_product expects.
inner(a::DifferentialFormValue{K,D,T,L},
      b::DifferentialFormValue{K,D,T,L}) where {K,D,T,L} =
  sum(a.data[i] * b.data[i] for i in 1:L)

# Direct indexing: bypass the generic MultiValue getindex which has an infinite
# recursion for rank-1 tensors (Int → CartesianIndex{1} → Int → ...).
@inline Base.getindex(v::DifferentialFormValue, i::Integer) = @inbounds v.data[i]

"""
    grad_to_2form(Jt::TensorValue{D,D})

Exterior derivative of a proxied 1-form `ω::VectorValue` from its gradient `Jt
= ∇ω`:

``(dω)_{a<b} = ∂ω_b/∂x^a − ∂ω_a/∂x^b = Jt[a,b] − Jt[b,a].``

Used to compute the exterior derivative of a vector proxied 1-form,, see `d_1form` in `Gridap.CellData`.
"""
function grad_to_2form(Jt::TensorValue{D,D,T,L}) where {D,T,L}
  cs = sorted_combinations(D, 2)
  DifferentialFormValue{2,D}(
    ntuple(n -> Jt[cs[n][1], cs[n][2]] - Jt[cs[n][2], cs[n][1]], binomial(D, 2)))
end
