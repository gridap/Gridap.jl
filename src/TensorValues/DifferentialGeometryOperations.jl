# Tensor operations involving proper tensor calculus and differential geometry operations.
#
# This includes pointwise algebra such as
# (∧, ι, ⋆, ♭/♯, koszul, …).

# ============================================================
# Generate-time expression builders
# ============================================================

# `sgn * expr`, with the sign folded in rather than emitted as a factor.
_signed(sgn::Int, expr) = sgn > 0 ? expr : :(-$expr)

# Sum of `terms`. `Expr(:call, :+)` with no term is `+()`, which errors, so an
# empty sum is emitted as an explicit `zero(T)`.
_sum_expr(terms, T) = isempty(terms) ? :(zero($T)) : Expr(:call, :+, terms...)

# Leibniz expansion of the determinant of the `K`×`K` matrix whose `(a,b)` entry
# is the expression `entry(a,b)`. `K ≥ 1`.
function _det_expr(entry, K::Int)
  terms = (_signed(levicivita(p), Expr(:call, :*, (entry(a, p[a]) for a in 1:K)...))
           for p in permutations(1:K))
  Expr(:call, :+, terms...)
end

# Determinant of the `K`×`K` submatrix of the second order tensor named by `A`,
# with rows `I` and columns `J`.
#
# The expansion has `K!` terms per call site. Past `K = 4` that makes compilation
# the dominant cost, so the emitted code calls `_minor_k` and loops instead.
function _minor_expr(A::Symbol, I, J)
  K = length(I)
  K == 0 && return :(one(eltype($A)))
  K > 4  && return :(_minor_k($A, $(Tuple(I)), $(Tuple(J))))
  _det_expr((a, b) -> :($A[$(I[a]), $(J[b])]), K)
end

# Determinant of the K×K submatrix of `A` with rows `I` and columns `J`.
function _minor_k(A, I::NTuple{K,Int}, J::NTuple{K,Int}) where K
  K == 0 && return one(eltype(A))
  det(SMatrix{K,K}(ntuple(n -> A[I[(n-1)%K+1], J[(n-1)÷K+1]], Val(K*K))))
end

# ============================================================
# Exterior (wedge) product ∧
# ============================================================

"""
    ∧(a::DifferentialFormValue{K1,D}, b::DifferentialFormValue{K2,D})

Pointwise exterior (wedge) product, a `DifferentialFormValue{K1+K2,D}`.
"""
@generated function ∧(a::DifferentialFormValue{K1,D,T1}, b::DifferentialFormValue{K2,D,T2}) where {K1,K2,D,T1,T2}
  K = K1 + K2
  L = binomial(D, K)
  T = Base.promote_op(*, T1, T2)

  iszero(L) && return :( zero(DifferentialFormValue{$K,$D,$T}) )

  # (a ∧ b)_I = Σ ε(I₁,I₂) a_{I₁} b_{I₂}, summed over the ways of splitting the
  # sorted I into a K1-subset I₁ and a K2-subset I₂. Overlapping splits have
  # ε = 0 and contribute nothing.
  terms = [Expr[] for _ in 1:L]
  for (n1, I1) in enumerate(sorted_combinations(D, K1)),
      (n2, I2) in enumerate(sorted_combinations(D, K2))
    s = sorting_sign(I1..., I2...)
    s == 0 && continue
    l = combination_index(sort([I1..., I2...]), D)
    push!(terms[l], _signed(s, :(indep_comp_getindex(a,$n1) * indep_comp_getindex(b,$n2))))
  end
  comps = [_sum_expr(t, T) for t in terms]

  :( DifferentialFormValue{$K,$D}(($(comps...),)) )
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
@generated function interior_product(v::VectorValue{D,Tv}, ω::DifferentialFormValue{K,D,Tw}) where {K,D,Tv,Tw}
  K >= 1 || return :(@unreachable "interior product requires K ≥ 1")

  Km1 = K - 1
  T = Base.promote_op(*, Tv, Tw)
  iszero(binomial(D, K)) && return :( zero(DifferentialFormValue{$Km1,$D,$T}) )

  L = binomial(D, Km1)

  # (ι_v ω)_{I∖I[j]} = Σ_j (-1)^{j-1} v^{I[j]} ω_I. Dropping one entry of the
  # sorted I leaves it sorted, so the output slot is its combination index.
  terms = [Expr[] for _ in 1:L]
  for (n, I) in enumerate(sorted_combinations(D, K)), j in 1:K
    l = combination_index(deleteat!(copy(I), j), D)
    push!(terms[l], _signed(iseven(j-1) ? 1 : -1,
                            :(v[$(I[j])] * indep_comp_getindex(ω,$n))))
  end
  comps = [_sum_expr(t, T) for t in terms]

  :( DifferentialFormValue{$Km1,$D}(($(comps...),)) )
end

const ι = interior_product

# ============================================================
# Hodge star ⋆ (flat Euclidean metric, g = Iᴅ, so √det g = 1)
#
# (⋆ω)_J = ∑_I  ε_{IJ}  ω_I
# where I ranges over K-combinations, J over (D-K)-combinations,
# and ε_{IJ} = levicivita of the permutation (I,J) → sorted (1…D).
# Only the complementary pair J = Iᶜ contributes, since ε_{IJ} vanishes as soon
# as I and J share an index.
# ============================================================

"""
    _hodge_star_signs(K, D) -> Vector{Int}

Levi-Civita signs of the flat Hodge star, entry `[n] = ε_{I Iᶜ}` for `I` the
n-th `K`-combination of `1:D`, so that with `L = binomial(D,K)`

    (⋆ω)_{L+1-n} = _hodge_star_signs(K,D)[n] * ω_n.

Only the complementary pair contributes, and its output slot is `L+1-n` because
complementation reverses the lexicographic order of the combinations. The sign
is `(-1)^(ΣI - K(K+1)/2)`, since the sorted `I` followed by its sorted
complement has `Σₐ (I[a] - a)` inversions.
"""
function _hodge_star_signs(K::Int, D::Int)

  # Helper: advance `I` to the next K-combination of 1:D in lexicographic
  # order, leaving the last combination unchanged.
  function _next_combination!(I::Vector{Int}, D::Int)
    K = length(I)
    j = K
    while j >= 1 && I[j] == D - K + j
      j -= 1
    end
    j == 0 && return I
    I[j] += 1
    for l in j+1:K
      I[l] = I[l-1] + 1
    end
    I
  end

  L = binomial(D, K)
  signs = Vector{Int}(undef, L)
  I = collect(1:K)              # first K-combination in lexicographic order
  offset = (K*(K+1)) ÷ 2
  for n in 1:L
    signs[n] = iseven(sum(I) - offset) ? 1 : -1
    _next_combination!(I, D)
  end
  signs
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
@generated function hodge_star(ω::DifferentialFormValue{K,D,T}) where {K,D,T}
  K <= D || return :(@unreachable $("hodge star requires K ≤ D, got K = $K and D = $D"))

  Kc = D - K
  L  = binomial(D, K)   # = binomial(D, Kc), the two are matched by complementation
  s  = _hodge_star_signs(K, D)
  comps = [_signed(s[L+1-m], :(indep_comp_getindex(ω,$(L+1-m)))) for m in 1:L]
  :( DifferentialFormValue{$Kc,$D}(($(comps...),)) )
end

const ⋆ = hodge_star

# ============================================================
# Generic functions extended elsewhere
#
# Methods for Symbolics.Num coefficients are provided by the GridapSymbolicsExt
# package extension.
# ============================================================

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

@generated function hodge_star(ω::DifferentialFormValue{K,D,Tw}, g_inv::SymTensorValue{D,Tg}, sqrt_det_g::Ts) where {K,D,Tw,Tg,Ts}
  K <= D || return :(@unreachable $("hodge star requires K ≤ D, got K = $K and D = $D"))

  T = Base.promote_op(*, Tw, Tg, Ts)
  Kc  = D - K
  L   = binomial(D, K)   # = binomial(D, Kc), the two are matched by complementation
  c_K = sorted_combinations(D, K)

  # Step 1: raise K indices of ω
  raised = [_sum_expr([:( $(_minor_expr(:g_inv, I, J)) * indep_comp_getindex(ω,$m) )
                       for (m, J) in enumerate(c_K)], T)
            for I in c_K]

  # Step 2: Levi-Civita contraction, which pairs each I with its complement only
  s = _hodge_star_signs(K, D)
  comps = [_signed(s[L+1-m], :(sqrt_det_g * ω_raised[$(L+1-m)])) for m in 1:L]

  quote
    ω_raised = ($(raised...),)
    DifferentialFormValue{$Kc,$D}(($(comps...),))
  end
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
@generated function flat(v::VectorValue{D,Tv}, g::SymTensorValue{D,Tg}) where {D,Tv,Tg}
  iszero(D) && return :( zero(DifferentialFormValue{1,0,$(Base.promote_op(*,Tg,Tv))}) )

  comps = [Expr(:call, :+, (:(g[$i,$j] * v[$j]) for j in 1:D)...) for i in 1:D]
  :( DifferentialFormValue{1,$D}(($(comps...),)) )
end

const ♭ = flat

"""
    sharp(ω::DifferentialFormValue{1,D}, g_inv::SymTensorValue{D})
    ♯(ω, g_inv)

Musical isomorphism `♯`: raise the index of the 1-form `ω` with the inverse
metric `g_inv`, giving a `VectorValue{D}`.
"""
@generated function sharp(ω::DifferentialFormValue{1,D,Tw}, g_inv::SymTensorValue{D,Tg}) where {D,Tw,Tg}
  iszero(D) && return :( zero(VectorValue{0,$(Base.promote_op(*,Tg,Tw))}) )

  comps = [Expr(:call, :+, (:(g_inv[$i,$j] * indep_comp_getindex(ω,$j)) for j in 1:D)...) for i in 1:D]
  :( VectorValue{$D}(($(comps...),)) )
end

const ♯ = sharp

"""
    apply_form(ω::DifferentialFormValue{K,D,T}, v₁, …, vₖ) → scalar

Evaluate the K-form ω on K tangent vectors v₁,…,vₖ ∈ ℝᴰ:
  ω(v₁,…,vₖ) = Σ_{|I|=K} ω_I · det([vₐ[I[b]]]_{a,b=1}^K)

Each vector must be a size `(D, )` tensor. Behaves like `getindex` for scalar `K=0`.
"""
@generated function apply_form(ω::DifferentialFormValue{K,D,T},
                               vs::Vararg{MultiValue{Tuple{D}},K}) where {K,D,T}
  K == 0 && return :(ω[1])

  # trivial Λᴷ for K > D
  iszero(binomial(D, K)) &&
    return Expr(:call, :*, :(zero($T)), (:(zero(eltype(vs[$a]))) for a in 1:K)...)

  terms = [:( indep_comp_getindex(ω,$n) * $(_det_expr((a,b) -> :(vs[$a][$(I[b])]), K)) )
           for (n, I) in enumerate(sorted_combinations(D, K))]
  Expr(:call, :+, terms...)
end

apply_form(ω::DifferentialFormValue{K,D}, vs...) where {K,D} =
  @unreachable "a K-form takes K vectors (subtyping `MultiValue`) of size (D,), got $(typeof.(vs)) for K = $K and D = $D"

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
@generated function pullback(ω::DifferentialFormValue{K,Dn,Tw}, J::TensorValue{Dn,Dm,Tj}) where {K,Dn,Dm,Tw,Tj}
  T = Base.promote_op(*, Tw, Tj)

  (iszero(binomial(Dn, K)) || iszero(binomial(Dm, K))) &&
    return :( zero(DifferentialFormValue{$K,$Dm,$T}) )

  c_I = sorted_combinations(Dn, K)   # K-combins in ambient space
  c_J = sorted_combinations(Dm, K)   # K-combins in chart space

  comps = [_sum_expr([:( $(_minor_expr(:J, I, Jidx)) * indep_comp_getindex(ω,$n) )
                      for (n, I) in enumerate(c_I)], T)
           for Jidx in c_J]

  :( DifferentialFormValue{$K,$Dm}(($(comps...),)) )
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
vol_coeff(ω::DifferentialFormValue{D,D}) where D = indep_comp_getindex(ω, 1)

# ============================================================
# Gradient support: outer(x, ω) for 1-forms.
#
# Gridap derives the gradient value type of a V-valued basis as
# typeof(outer(zero(x), zero(V))) (Fields.gradient_type). Defining outer of a
# point with a 1-form value as the TensorValue{D,D} with entries
# [i,j] = x_i·ω_j makes ∇(basis) follow Gridap's convention ∇u[i,j] = ∂u_j/∂xⁱ
# for the vector proxy of the 1-form.
# ============================================================

@generated function outer(a::VectorValue{D,Ta}, b::DifferentialFormValue{1,D,Tb}) where {D,Ta,Tb}
  iszero(D) && return :( zero(TensorValue{0,0,$(Base.promote_op(*,Ta,Tb))}) )

  comps = [:( a[$i] * indep_comp_getindex(b,$j) ) for j in 1:D for i in 1:D]
  :( TensorValue{$D,$D}(($(comps...),)) )
end

# ============================================================
# Inner product of K-forms
#
# (α|β) = (1/K!) α_{i₁…i_K} g^{i₁j₁} ⋯ g^{i_Kj_K} β_{j₁…j_K}
#       = Σ_{|I|=|J|=K} α_I β_J det([g⁻¹[I[a],J[b]]]_{a,b})
#
# the second sum running over sorted multi-indices, which reduces to
# Σ_I α_I β_I for a flat metric. This is the inner product the Hodge star is
# built on,
#   α ∧ ⋆_g β = (α|β) μ = (α|β) √det(g) dx¹ ∧ … ∧ dx^D.
#
# `inner` is the tensor inner product, such that α⊙β = K!(α|β) for a flat metric.
# ============================================================

"""
    form_inner(ω::DifferentialFormValue{K,D}, η::DifferentialFormValue{K,D})
    form_inner(ω::DifferentialFormValue{K,D}, η::DifferentialFormValue{K,D}, g_inv::SymTensorValue{D})
    ω ⨟ η

Inner product ``(ω|η)`` of two K-forms, the scalar such that

``ω ∧ ⋆_g η = √det(g) (ω|η) dx¹∧…∧dx^D``.

The two-argument form uses the flat Euclidean metric. The three-argument form
takes the inverse metric tensor.

This differs from the inner product of general tensors, `ω ⊙ η == factorial(K) (ω ⨟ μ)`.
"""
@generated function form_inner(a::DifferentialFormValue{K,D,Ta},
                               b::DifferentialFormValue{K,D,Tb}) where {K,D,Ta,Tb}
  T = Base.promote_op(*,Ta,Tb)
  L = num_indep_components(a)
  _sum_expr([:(indep_comp_getindex(a,$i) * indep_comp_getindex(b,$i)) for i in 1:L], T)
end

@generated function form_inner(a::DifferentialFormValue{K,D,Ta}, b::DifferentialFormValue{K,D,Tb},
                               g_inv::SymTensorValue{D,Tg}) where {K,D,Ta,Tb,Tg}
  T = Base.promote_op(*, Ta, Tb, Tg)
  c_K = sorted_combinations(D, K)
  # raising β_I to β^I as `hodge_star` does
  terms = [:( indep_comp_getindex(a,$m) *
              $(_sum_expr([:( $(_minor_expr(:g_inv, I, J)) * indep_comp_getindex(b,$n) )
                           for (n, J) in enumerate(c_K)], T)) )
           for (m, I) in enumerate(c_K)]
  # (α|β) = Σ_I α_I β^I
  _sum_expr(terms, T)
end

const ⨟ = form_inner

"""
    grad_to_2form(Jt::TensorValue{D,D})

Exterior derivative of a proxied 1-form `ω::VectorValue` from its gradient `Jt
= ∇ω`:

``(dω)_{a<b} = ∂ω_b/∂x^a − ∂ω_a/∂x^b = Jt[a,b] − Jt[b,a].``

Used to compute the exterior derivative of a vector proxied 1-form,, see `d_1form` in `Gridap.CellData`.
"""
@generated function grad_to_2form(Jt::TensorValue{D,D,T}) where {D,T}
  iszero(binomial(D, 2)) && return :( zero(DifferentialFormValue{2,$D,$(Base.promote_op(-,T,T))}) )

  comps = [:( Jt[$a,$b] - Jt[$b,$a] ) for (a, b) in sorted_combinations(D, 2)]
  :( DifferentialFormValue{2,$D}(($(comps...),)) )
end


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
to_1form(v::VectorValue{D,T}) where {D,T} = DifferentialFormValue{1,D,T}(Tuple(v))

"""
    from_1form(ω::DifferentialFormValue{1})

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
to_Kform(v::VectorValue{L,T}, ::Val{K}, ::Val{D}) where {L,T,K,D} = DifferentialFormValue{K,D,T}(Tuple(v))

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

