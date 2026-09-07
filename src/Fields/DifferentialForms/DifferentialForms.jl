# DifferentialForms.jl
#
# Lazy Field-level differential forms, following the Gridap Field pattern:
#   - DifferentialForm{K,D,L,Data}  — typed data, analogous to FieldGradient{O,F}
#   - ExteriorDerivativeForm{K,D,F} — lazy wrapper, stores only form::F
#   - CodifferentialForm{K,D,F}     — lazy wrapper, stores only form::F
# All computation deferred to return_cache / evaluate!.
#
# CRITICAL: never use lambda-wrapped OperationField (e.g., Operation(x->-x)(f))
# in component fields.  FieldGradient cannot differentiate those.
# Use the Field arithmetic that Gridap defines:
#   -(f::Field)           → Operation(-)(f)  — gradient defined
#   s * f                 → Operation(*)(ConstantField(s), f)
#   Operation(+)(f,g)     — gradient defined

using Combinatorics: levicivita

import Gridap.TensorValues: ∧
import Gridap.TensorValues: koszul
import Gridap.TensorValues: pullback

# ============================================================
# DifferentialForm{K,D,L,Data}
# ============================================================

"""
    DifferentialForm{K,D,L,Data} <: Field

Typed K-form in D dimensions.  `Data = typeof(data)` captures the full tuple
type of the L = binomial(D,K) component Fields so that downstream wrappers
(ExteriorDerivativeForm, CodifferentialForm) can use `map(∇, form.data)` with
full type inference — giving zero-allocation, type-stable evaluation.
"""
struct DifferentialForm{K,D,L,Data} <: Field
  data :: Data
  function DifferentialForm{K,D}(data) where {K,D}
    L = length(data)
    @assert L == binomial(D,K) "wrong number of component fields"
    new{K,D,L,typeof(data)}(data)
  end
end

testargs(a::DifferentialForm, x::Point) = testargs.(a.data, Ref(x))

function return_value(a::DifferentialForm{K,D}, x::Point) where {K,D}
  DifferentialFormValue{K,D}(map(f -> return_value(f, x), a.data))
end

function return_cache(a::DifferentialForm, x::Point)
  map(f -> return_cache(f, x), a.data)   # typed tuple of caches
end

function evaluate!(cache, a::DifferentialForm{K,D}, x::Point) where {K,D}
  DifferentialFormValue{K,D}(map((c,f) -> evaluate!(c,f,x), cache, a.data))
end


# ============================================================
# ExteriorDerivativeForm{K,D,F} — lazy exterior derivative
#
# Analogous to FieldGradient{O,F}: stores only the wrapped form,
# defers gradient computation to return_cache / evaluate!.
#
# return_cache uses map(∇, form.data) which is type-stable when Data is typed.
# evaluate! uses map over the typed tuple of (cache, grad_field) pairs
# — zero allocation for bits-type DifferentialFormValue.
# ============================================================

"""
    ExteriorDerivativeForm{K,D,F} <: Field

Lazy exterior derivative of a `DifferentialForm{K,D,...}`.  Stores only
`form::F`; gradient fields are derived in `return_cache` via `map(∇, form.data)`.
"""
struct ExteriorDerivativeForm{K,D,F} <: Field
  form :: F
  function ExteriorDerivativeForm(form::DifferentialForm{K,D}) where {K,D}
    new{K,D,typeof(form)}(form)
  end
end

exterior_derivative(f::DifferentialForm) = ExteriorDerivativeForm(f)

function return_cache(dω::ExteriorDerivativeForm{K,D,F}, x::Point) where {K,D,F}
  # map(∇, form.data) is type-stable when Data = Tuple{F1,F2,...}
  grad_fields = map(∇, dω.form.data)
  caches      = map(g -> return_cache(g, x), grad_fields)
  (grad_fields, caches)
end

function return_value(dω::ExteriorDerivativeForm{K,D,F}, x::Point) where {K,D,F}
  evaluate(dω, x)
end

function evaluate!(cache, dω::ExteriorDerivativeForm{K,D,F}, x::Point) where {K,D,F}
  grad_fields, caches = cache
  L = binomial(D, K)
  # Evaluate all gradients via map — type-stable, no Vector allocation
  gvs = map((c,g) -> evaluate!(c, g, x), caches, grad_fields)
  # Accumulate using an explicit loop (map-reduce over tuple is type-stable)
  gv1  = gvs[1]
  dfi1 = DifferentialFormValue{1,D}(gv1.data)
  kbi1 = DifferentialFormValue{K,D}(ntuple(j -> j == 1 ? 1.0 : 0.0, Val(L)))
  acc  = dfi1 ∧ kbi1
  for i in 2:L
    gv  = gvs[i]
    dfi = DifferentialFormValue{1,D}(gv.data)
    kbi = DifferentialFormValue{K,D}(ntuple(j -> j == i ? 1.0 : 0.0, Val(L)))
    acc = acc + (dfi ∧ kbi)
  end
  acc
end

# ============================================================
# hodge_star_form: Field-level flat Hodge star
#
# Builds DifferentialForm{D-K,D} with component fields that are ±1 linear
# combinations of the original components.
# IMPORTANT: uses -(f::Field) = Operation(-)(f) for negation.
# This has a defined gradient rule (gradient of sum/difference),
# unlike the lambda Operation(x->-x)(f) which hits @abstractmethod.
# ============================================================

"""
    hodge_star_form(ω::DifferentialForm{K,D})

Flat (Euclidean) Hodge star: a `DifferentialForm{D-K,D}` whose component fields
are ±1 linear combinations of the components of `ω`.

It is differentiable, unlike the pointwise `Operation(hodge_star)`.
"""
function hodge_star_form(ω::DifferentialForm{K,D}) where {K,D}
  Kc = D - K
  L  = binomial(D, K)
  Lc = binomial(D, Kc)

  c_in  = sorted_combinations(D, K)
  c_out = sorted_combinations(D, Kc)

  component_fields = ntuple(Lc) do m
    J = c_out[m]
    accum = nothing
    for (n, I) in enumerate(c_in)
      perm = [I..., J...]
      length(unique(perm)) == D || continue
      sgn  = levicivita(sortperm(perm))
      # KEY: use -(f) not Operation(x->-x)(f) — the former is differentiable
      term = sgn == 1 ? ω.data[n] : -(ω.data[n])
      accum = isnothing(accum) ? term : Operation(+)(accum, term)
    end
    accum
  end

  DifferentialForm{Kc,D}(component_fields)
end

# ExteriorDerivativeForm evaluates to DifferentialFormValue{K+1,D}.
# hodge_star_form on it: extract component fields via _component_field.
# NOTE: the resulting component fields use Operation(x->x.data[i]) lambdas.
# These are NOT differentiable by FieldGradient.
# This method is safe for evaluation only — do NOT apply exterior_derivative
# to its result.
_component_field(f::Field, i::Int) = Operation(x -> x.data[i])(f)

function hodge_star_form(ω::ExteriorDerivativeForm{K,D,F}) where {K,D,F}
  Lp1   = binomial(D, K+1)
  comps = ntuple(i -> _component_field(ω, i), Lp1)
  hodge_star_form(DifferentialForm{K+1,D}(comps))
end

# ============================================================
# CodifferentialForm{K,D,F} — lazy codifferential
#   δω = (-1)^{D(K-1)+1} ⋆ d ⋆ ω
#
# Stores only form::F.  return_cache derives grad_fields and Hodge-star
# matrix; evaluate! is purely algebraic — zero allocation.
# ============================================================

"""
    CodifferentialForm{K,D,F} <: Field

Lazy codifferential of `DifferentialForm{K,D,...}`.  Stores only `form::F`;
all helper data (Hodge star matrix, gradient caches) are derived in
`return_cache` and reused across point evaluations.
"""
struct CodifferentialForm{K,D,F} <: Field
  form :: F
  function CodifferentialForm(form::DifferentialForm{K,D}) where {K,D}
    @assert K >= 1 "codifferential is zero on 0-forms"
    new{K,D,typeof(form)}(form)
  end
end

codifferential(ω::DifferentialForm{K,D}) where {K,D} = CodifferentialForm(ω)

function return_cache(δω::CodifferentialForm{K,D,F}, x::Point) where {K,D,F}
  L   = binomial(D, K)
  Lc  = binomial(D, D-K)

  # Hodge star sign matrix (constant — computed once per cache setup)
  hs_mat = TensorValues._hodge_star_matrix(K, D)
  sgn = iseven(D*(K-1) + 1) ? 1 : -1

  # Gradient fields and their caches — typed via map(∇, form.data)
  grad_fields = map(∇, δω.form.data)
  grad_caches = map(g -> return_cache(g, x), grad_fields)

  (sgn, hs_mat, grad_fields, grad_caches)
end

function return_value(δω::CodifferentialForm{K,D,F}, x::Point) where {K,D,F}
  evaluate(δω, x)
end

function evaluate!(cache, δω::CodifferentialForm{K,D,F}, x::Point) where {K,D,F}
  sgn, hs_mat, grad_fields, grad_caches = cache
  L   = binomial(D, K)
  Lc  = binomial(D, D-K)

  # Evaluate all gradients — map over typed tuple, no Vector allocation
  gvs = map((c,g) -> evaluate!(c, g, x), grad_caches, grad_fields)

  # ∇(⋆ω)_m = Σ_n hs_mat[m,n] · gvs[n]   (VectorValue{D})
  star_grads = ntuple(Val(Lc)) do m
    s = hs_mat[m, 1] * gvs[1]
    for n in 2:L
      s = s + hs_mat[m, n] * gvs[n]
    end
    s
  end

  # d(⋆ω): same accumulation as ExteriorDerivativeForm but using star_grads
  dfi1 = DifferentialFormValue{1,D}(star_grads[1].data)
  kbm1 = DifferentialFormValue{D-K,D}(ntuple(j -> j == 1 ? 1.0 : 0.0, Val(Lc)))
  acc  = dfi1 ∧ kbm1
  for m in 2:Lc
    dfi = DifferentialFormValue{1,D}(star_grads[m].data)
    kbm = DifferentialFormValue{D-K,D}(ntuple(j -> j == m ? 1.0 : 0.0, Val(Lc)))
    acc = acc + (dfi ∧ kbm)
  end

  sgn * hodge_star(acc)   # DifferentialFormValue{K-1,D}
end


# ============================================================
# KoszulForm{K,D,F} — lazy Koszul contraction  κ_x(ω)
#
# (κω)(x) = ι_x(ω(x))  — interior product with the evaluation point.
# Since Point{D,T} = VectorValue{D,T}, x can be passed
# directly to interior_product.
#
# Homotopy identity:  dκ + κd = (r+K) id  on Pᵣ Λᴷ.
# This is the algebraic kernel of trimmed polynomial form spaces.
# ============================================================

"""
    KoszulForm{K,D,F} <: Field

Lazy Koszul contraction of a `DifferentialForm{K,D,...}`.
Stores only `form::F`; evaluates `ι_x(form(x))` at each point x.
"""
struct KoszulForm{K,D,F} <: Field
  form :: F
  function KoszulForm(form::DifferentialForm{K,D}) where {K,D}
    @assert K >= 1 "Koszul operator requires K ≥ 1"
    new{K,D,typeof(form)}(form)
  end
end

"""
    koszul(ω::DifferentialForm)

Koszul differential of `ω`, giving a (K-1)-form valued `KoszulForm`.
"""
koszul(form::DifferentialForm) = KoszulForm(form)

function return_cache(κω::KoszulForm{K,D,F}, x::Point) where {K,D,F}
  return_cache(κω.form, x)
end

function return_value(κω::KoszulForm{K,D,F}, x::Point) where {K,D,F}
  evaluate(κω, x)
end

function evaluate!(cache, κω::KoszulForm{K,D,F}, x::Point) where {K,D,F}
  ω_x = evaluate!(cache, κω.form, x)
  interior_product(x, ω_x)   # Point{D,T} = VectorValue{D,T}
end

# ============================================================
# PullbackForm{K,Dm,Dn} — Field-level pullback  φ*ω
# ============================================================

"""
    PullbackForm{K,Dm,Dn} <: Field

Lazy pullback `φ*ω` of a `DifferentialForm{K,Dn}` under a map field
`φ : ℝᴰᵐ → ℝᴰⁿ`, see [`pullback(φ::Field, ω::DifferentialForm, ::Val)`](@ref pullback).
"""
struct PullbackForm{K,Dm,Dn} <: Field
  map_field :: Any
  jac_field :: Any
  form      :: Any   # DifferentialForm{K,Dn,...}
end

function pullback(φ::Field, ω::DifferentialForm{K,Dn}, ::Val{Dm}) where {K,Dn,Dm}
  PullbackForm{K,Dm,Dn}(φ, ∇(φ), ω)
end

function return_cache(f::PullbackForm{K,Dm,Dn}, x::Point{Dm}) where {K,Dm,Dn}
  mc = return_cache(f.map_field, x)
  y0 = evaluate!(mc, f.map_field, x)
  fc = return_cache(f.form, y0)
  jc = return_cache(f.jac_field, x)
  (mc, fc, jc)
end

function return_value(f::PullbackForm{K,Dm,Dn}, x::Point{Dm}) where {K,Dm,Dn}
  evaluate(f, x)
end

function evaluate!(cache, f::PullbackForm{K,Dm,Dn}, x::Point{Dm}) where {K,Dm,Dn}
  mc, fc, jc = cache
  y   = evaluate!(mc, f.map_field, x)
  J   = evaluate!(jc, f.jac_field, x)
  ω_y = evaluate!(fc, f.form, y)
  pullback(ω_y, J)
end

# ============================================================
# Vector-of-points evaluation, shared by all lazy form Fields
# ============================================================

const _LazyFormField = Union{DifferentialForm,ExteriorDerivativeForm,
                             CodifferentialForm,KoszulForm,PullbackForm}

function return_cache(f::_LazyFormField, x::AbstractVector{<:Point})
  xi = testitem(x)
  cf = return_cache(f, xi)
  vi = evaluate!(cf, f, xi)
  r  = CachedArray(Vector{typeof(vi)}(undef, length(x)))
  r, cf
end

function evaluate!(cache, f::_LazyFormField, x::AbstractVector{<:Point})
  r, cf = cache
  np = length(x)
  setsize!(r, (np,))
  for i in 1:np
    @inbounds r.array[i] = evaluate!(cf, f, x[i])
  end
  r.array
end
