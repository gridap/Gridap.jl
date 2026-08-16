# DifferentialFormCellFields.jl
#
# CellField support for differential forms.
#
# Design:
#  - DifferentialFormCellField{K,D,DS} <: CellField carries the form degree K
#    and ambient dimension D as type parameters.  The underlying cell array
#    stores any <:Field that evaluates to DifferentialFormValue{K,D}.
#
#  - Form-level ops (exterior_derivative, codifferential, hodge_star_form) use
#    lazy_map(Broadcasting(op), get_data(a)) and return a new
#    DifferentialFormCellField with the updated degree.  All computation is
#    deferred; the LazyArray reuses a single pre-allocated cache per cell.
#
#  - Value-level ops (∧, hodge_star, flat, sharp, ι) use Operation(op)(args...)
#    which creates an OperationCellField.  At evaluation time
#    BroadcastingFieldOpMap evaluates op(aᵢ, bᵢ) in a preallocated buffer —
#    zero allocation in the inner loop for bits-type DifferentialFormValue.

import Gridap.TensorValues: ∧
import Gridap.TensorValues: hodge_star
import Gridap.TensorValues: flat
import Gridap.TensorValues: sharp
import Gridap.TensorValues: interior_product
import Gridap.TensorValues: exterior_derivative
import Gridap.TensorValues: codifferential
import Gridap.TensorValues: koszul
import Gridap.TensorValues: vol_coeff
import Gridap.TensorValues: to_1form
import Gridap.TensorValues: from_1form
import Gridap.Fields: hodge_star_form

# ============================================================
# DifferentialFormCellField
# ============================================================

"""
    DifferentialFormCellField{K,D,DS} <: CellField

Cell field of K-forms in D dimensions.
- `cell_field`: array whose i-th entry is a `Field` that evaluates to
  `DifferentialFormValue{K,D}` at a given point.
- `trian`: underlying `Triangulation`
- `domain_style`: `ReferenceDomain()` or `PhysicalDomain()`
"""
struct DifferentialFormCellField{K,D,DS} <: CellField
  cell_field   :: AbstractArray
  trian        :: Triangulation
  domain_style :: DS

  function DifferentialFormCellField{K,D}(
      cell_field   :: AbstractArray,
      trian        :: Triangulation,
      domain_style :: DomainStyle) where {K,D}
    DS = typeof(domain_style)
    new{K,D,DS}(Fields.MemoArray(cell_field), trian, domain_style)
  end
end

# ── CellDatum / CellField interface ──────────────────────────────────────────

get_data(f::DifferentialFormCellField)          = f.cell_field
get_triangulation(f::DifferentialFormCellField) = f.trian
DomainStyle(::Type{DifferentialFormCellField{K,D,DS}}) where {K,D,DS} = DS()

# similar_cell_field preserves K and D through change_domain etc.
function similar_cell_field(
    f  :: DifferentialFormCellField{K,D},
    data, trian, ds) where {K,D}
  DifferentialFormCellField{K,D}(data, trian, ds)
end

# ── Convenience constructors ──────────────────────────────────────────────────

"""
    DifferentialFormCellField{K,D}(form::DifferentialForm{K,D}, trian, domain_style)

Uniform field: every cell shares the same `DifferentialForm` object.
"""
function DifferentialFormCellField{K,D}(
    form  :: DifferentialForm{K,D},
    trian :: Triangulation,
    domain_style :: DomainStyle) where {K,D}
  cell_field = Fill(form, num_cells(trian))
  DifferentialFormCellField{K,D}(cell_field, trian, domain_style)
end

# Note: there is intentionally no Function constructor.
# exterior_derivative (and codifferential, etc.) need the component fields
# of a DifferentialForm to be accessible individually for gradient computation.
# A GenericField black-box function does not expose those components.
# Construct from DifferentialForm (or an AbstractArray of them) instead.

# ── Form-level operators (lazy_map, zero allocation per cell) ─────────────────

"""
    exterior_derivative(a::DifferentialFormCellField{K,D})

Returns a `DifferentialFormCellField{K+1,D}` whose cell data are
`ExteriorDerivativeForm{K,D}` objects (pre-built at construction,
gradient fields pre-allocated — evaluate! is allocation-free).
"""
function exterior_derivative(a::DifferentialFormCellField{K,D}) where {K,D}
  cell_da = lazy_map(Broadcasting(exterior_derivative), get_data(a))
  DifferentialFormCellField{K+1,D}(cell_da, get_triangulation(a), DomainStyle(a))
end

"""
    codifferential(a::DifferentialFormCellField{K,D})

Returns a `DifferentialFormCellField{K-1,D}`.  Each cell's
`CodifferentialForm` precomputes the Hodge star matrix and gradient caches
at construction; `evaluate!` is allocation-free.
"""
function codifferential(a::DifferentialFormCellField{K,D}) where {K,D}
  @assert K >= 1 "codifferential is zero on 0-forms"
  cell_δa = lazy_map(Broadcasting(codifferential), get_data(a))
  DifferentialFormCellField{K-1,D}(cell_δa, get_triangulation(a), DomainStyle(a))
end

"""
    hodge_star_form(a::DifferentialFormCellField{K,D})

Returns a `DifferentialFormCellField{D-K,D}` by applying the flat
Euclidean `hodge_star_form` to each component field (constant ±1
linear combinations — no spatial dependence in the coefficients).
Valid when the cell data is `DifferentialForm` objects.
"""
function hodge_star_form(a::DifferentialFormCellField{K,D}) where {K,D}
  cell_star = lazy_map(Broadcasting(hodge_star_form), get_data(a))
  DifferentialFormCellField{D-K,D}(cell_star, get_triangulation(a), DomainStyle(a))
end

# ── Value-level operators via Operation (BroadcastingFieldOpMap, zero alloc) ──
# These defer to the OperationCellField machinery.  At evaluation time,
# BroadcastingFieldOpMap applies op(aᵢ, bᵢ) into a preallocated output array.

"""
    (∧)(a::DifferentialFormCellField, b::DifferentialFormCellField)

Pointwise exterior product.  Returns an `OperationCellField`.
"""
∧(a::DifferentialFormCellField, b::DifferentialFormCellField) = Operation(∧)(a, b)

# Accept mixing with plain CellField (e.g., scalar CellField ∧ form)
∧(a::CellField, b::DifferentialFormCellField) = Operation(∧)(a, b)
∧(a::DifferentialFormCellField, b::CellField) = Operation(∧)(a, b)

# Fallback: two generic CellFields (e.g., both OperationCellField from to_1form/d_0form)
∧(a::CellField, b::CellField) = Operation(∧)(a, b)

"""
    hodge_star(a::DifferentialFormCellField{K,D})

Pointwise flat Hodge star (g = I_D).  Returns an `OperationCellField`.
"""
hodge_star(a::DifferentialFormCellField) = Operation(hodge_star)(a)

# Fallback: generic CellField (e.g., OperationCellField from to_1form/d_0form)
hodge_star(a::CellField) = Operation(hodge_star)(a)

"""
    flat(v::CellField, g::CellField)

Lower a vector field `v` using metric `g` (a `SymTensorValue`-valued CellField).
"""
flat(v::CellField, g::CellField) = Operation(flat)(v, g)

"""
    sharp(ω::DifferentialFormCellField, g_inv::CellField)

Raise a 1-form `ω` using inverse metric `g_inv`.
"""
sharp(ω::DifferentialFormCellField, g_inv::CellField) = Operation(sharp)(ω, g_inv)

"""
    interior_product(v::CellField, ω::DifferentialFormCellField)
"""
interior_product(v::CellField, ω::DifferentialFormCellField) = Operation(interior_product)(v, ω)

# ── Arithmetic ────────────────────────────────────────────────────────────────

Base.:+(a::DifferentialFormCellField, b::DifferentialFormCellField) = Operation(+)(a, b)
Base.:-(a::DifferentialFormCellField, b::DifferentialFormCellField) = Operation(-)(a, b)
Base.:-(a::DifferentialFormCellField) = Operation(-)(a)

function Base.:*(s::Union{Real,Complex}, a::DifferentialFormCellField)
  Operation(Base.:*)(CellField(s, get_triangulation(a), DomainStyle(a)), a)
end

Base.:*(a::DifferentialFormCellField, s::Union{Real,Complex}) = s * a

# ── Discrete exterior derivative for FEFunctions ──────────────────────────────
#
# When u comes from a Gridap FESpace (not a DifferentialForm), its component
# fields are not directly accessible, so exterior_derivative(DifferentialForm)
# cannot be applied.  These functions compute d(u) from the gradient tensor
# using algebraic identities — no component access required.
#
#   d_0form(u): u is scalar-valued (Lagrange FEFunction, 0-form)
#     du = to_1form(∇u) = Σᵢ (∂u/∂xⁱ) dxⁱ
#
#   d_1form(u): u is VectorValue-valued (Nédélec FEFunction, 1-form)
#     du = jac_to_2form(∇u),  (dω)_{a<b} = ∂ω_b/∂x^a − ∂ω_a/∂x^b

"""
    d_0form(u::CellField)

Discrete exterior derivative of a scalar (0-form) CellField.
Returns a `DifferentialFormValue{1,D}`-valued `OperationCellField`.

Equivalent to `to_1form(∇(u))`.
"""
d_0form(u::CellField) = to_1form(∇(u))

"""
    d_1form(u::CellField)

Discrete exterior derivative of a VectorValue (1-form) CellField
(e.g., from a Nédélec FESpace).
Returns a `DifferentialFormValue{2,D}`-valued `OperationCellField`.

Computed as `Operation(jac_to_2form)(∇(u))`.
"""
d_1form(u::CellField) = Operation(jac_to_2form)(∇(u))

# ── Koszul operator (lazy_map, zero allocation per cell) ──────────────────────

"""
    koszul(a::DifferentialFormCellField{K,D})

Returns a `DifferentialFormCellField{K-1,D}` by applying `KoszulForm` to
each cell's `DifferentialForm`.  Evaluation contracts `form(x)` with the
point x at each quadrature node — no allocation in the inner loop.
"""
function koszul(a::DifferentialFormCellField{K,D}) where {K,D}
  @assert K >= 1 "Koszul operator requires K ≥ 1"
  cell_κa = lazy_map(Broadcasting(koszul), get_data(a))
  DifferentialFormCellField{K-1,D}(cell_κa, get_triangulation(a), DomainStyle(a))
end

# ── Volume form scalar extraction ──────────────────────────────────────────────

"""
    vol_coeff(a::DifferentialFormCellField{D,D})

Extracts the single scalar coefficient of a top-degree D-form CellField.
Returns an `OperationCellField` (scalar-valued) suitable for use in
`∫(vol_coeff(ω ∧ ⋆η)) * dΩ` with Gridap's standard measure.
"""
vol_coeff(a::DifferentialFormCellField{D,D}) where {D} = Operation(vol_coeff)(a)
vol_coeff(a::CellField) = Operation(vol_coeff)(a)  # OperationCellField from ω ∧ ⋆η

# ── VectorValue ↔ DifferentialFormValue{1,D} bridge ───────────────────────────

"""
    to_1form(v::CellField)

Pointwise conversion of a VectorValue{D}-valued CellField (e.g., a Nédélec
FEFunction) to a DifferentialFormValue{1,D}-valued OperationCellField.
The result can be used with the value-level operators (∧, ⋆, ι).
Note: for `exterior_derivative`, wrap as `DifferentialFormCellField{1,D}`
with accessible component fields instead.
"""
to_1form(v::CellField) = Operation(to_1form)(v)

"""
    from_1form(ω::DifferentialFormCellField{1,D})

Pointwise conversion of a DifferentialFormValue{1,D}-valued CellField to a
VectorValue{D}-valued OperationCellField.
"""
from_1form(ω::DifferentialFormCellField{1,D}) where {D} = Operation(from_1form)(ω)
