# DifferentialFormCellFields.jl
#
# CellField-level support for differential forms.
#
# Design:
#  - There is deliberately no dedicated `CellField` subtype for forms.  The form
#    degree K and ambient dimension D are carried by the *cell data* — the
#    `DifferentialForm{K,D}` fields and the `DifferentialFormValue{K,D}` values
#    they produce — which enforce those parameters at construction.  A CellField
#    wrapper repeating K and D could only duplicate that information without
#    verifying it, so form cell fields are plain `GenericCellField`s.
#
#  - Form-level ops (exterior_derivative, codifferential, hodge_star_form,
#    koszul) use lazy_map(Broadcasting(op), get_data(a)); the returned CellField
#    is built with similar_cell_field.  All computation is deferred; the
#    LazyArray reuses a single pre-allocated cache per cell.  The type check that
#    the cell data really are forms is done by `_form_op` below, and the degree
#    invariants (K ≥ 1 etc.) are asserted by the Field-level constructors.
#
#  - Value-level ops (∧, hodge_star, flat, sharp, ι) use Operation(op)(args...)
#    which creates an OperationCellField.  At evaluation time
#    BroadcastingFieldOpMap evaluates op(aᵢ, bᵢ) in a preallocated buffer —
#    zero allocation in the inner loop for bits-type DifferentialFormValue.
#
#  - Arithmetic (+, -, scalar *) is inherited from the generic CellField
#    operators in CellFields.jl; no form-specific methods are needed.

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
using Gridap.Fields: DifferentialForm

# ============================================================
# Construction
# ============================================================

"""
    form_cell_field(form::DifferentialForm{K,D}, trian, domain_style)

Cell field of K-forms in D dimensions in which every cell shares the same
`DifferentialForm` object.  Returns a plain `GenericCellField`: the form degree
is carried by the cell data, not by the CellField type.

To build a cell field with per-cell forms, pass an array of `DifferentialForm`s
to `GenericCellField` directly.

Note: there is intentionally no `Function` constructor.  `exterior_derivative`
(and `codifferential`, etc.) need the component fields of a `DifferentialForm`
to be accessible individually for gradient computation; a `GenericField`
black-box function does not expose those components.
"""
function form_cell_field(
    form  :: DifferentialForm,
    trian :: Triangulation,
    domain_style :: DomainStyle)
  GenericCellField(Fill(form, num_cells(trian)), trian, domain_style)
end

# ============================================================
# Form-level operators (lazy_map, zero allocation per cell)
# ============================================================

# Guard replacing the former `DifferentialFormCellField` type barrier: verify
# that the cell data are form fields and report it clearly, rather than letting
# a bare MethodError escape from inside lazy_map.
function _form_op(op, a::CellField, opname::String)
  data = get_data(a)
  f = testitem(data)
  hasmethod(op, Tuple{typeof(f)}) || throw(ArgumentError(
    "$opname requires a CellField whose cell data are differential forms; got " *
    "cell data of type $(typeof(f)). Build it from `DifferentialForm` objects " *
    "(e.g. via `form_cell_field`)."))
  similar_cell_field(a, lazy_map(Broadcasting(op), data),
                     get_triangulation(a), DomainStyle(a))
end

"""
    exterior_derivative(a::CellField)

Exterior derivative of a K-form cell field, giving a (K+1)-form cell field whose
cell data are `ExteriorDerivativeForm` objects (gradient fields pre-allocated at
construction — `evaluate!` is allocation-free).
"""
exterior_derivative(a::CellField) = _form_op(exterior_derivative, a, "exterior_derivative")

"""
    codifferential(a::CellField)

Codifferential of a K-form cell field, giving a (K-1)-form cell field.  Each
cell's `CodifferentialForm` precomputes the Hodge star matrix and gradient
caches at construction; `evaluate!` is allocation-free.  Asserts K ≥ 1.
"""
codifferential(a::CellField) = _form_op(codifferential, a, "codifferential")

"""
    hodge_star_form(a::CellField)

Flat Euclidean Hodge star applied to each component field of a K-form cell
field, giving a (D-K)-form cell field (constant ±1 linear combinations — no
spatial dependence in the coefficients).
"""
hodge_star_form(a::CellField) = _form_op(hodge_star_form, a, "hodge_star_form")

"""
    koszul(a::CellField)

Koszul contraction of a K-form cell field, giving a (K-1)-form cell field.
Evaluation contracts `form(x)` with the point x at each quadrature node — no
allocation in the inner loop.  Asserts K ≥ 1.
"""
koszul(a::CellField) = _form_op(koszul, a, "koszul")

# ============================================================
# Value-level operators via Operation (BroadcastingFieldOpMap, zero alloc)
# ============================================================
# These defer to the OperationCellField machinery.  At evaluation time,
# BroadcastingFieldOpMap applies op(aᵢ, bᵢ) into a preallocated output array.

"""
    (∧)(a::CellField, b::CellField)

Pointwise exterior product.  Returns an `OperationCellField`.
"""
∧(a::CellField, b::CellField) = Operation(∧)(a, b)

"""
    hodge_star(a::CellField)

Pointwise flat Hodge star (g = I_D).  Returns an `OperationCellField`.
"""
hodge_star(a::CellField) = Operation(hodge_star)(a)

"""
    flat(v::CellField, g::CellField)

Lower a vector field `v` using metric `g` (a `SymTensorValue`-valued CellField).
"""
flat(v::CellField, g::CellField) = Operation(flat)(v, g)

"""
    sharp(ω::CellField, g_inv::CellField)

Raise a 1-form `ω` using inverse metric `g_inv`.
"""
sharp(ω::CellField, g_inv::CellField) = Operation(sharp)(ω, g_inv)

"""
    interior_product(v::CellField, ω::CellField)
"""
interior_product(v::CellField, ω::CellField) = Operation(interior_product)(v, ω)

# ============================================================
# Discrete exterior derivative for FEFunctions
# ============================================================
#
# When u comes from a Gridap FESpace (not a DifferentialForm), its component
# fields are not directly accessible, so exterior_derivative cannot be applied.
# These functions compute d(u) from the gradient tensor using algebraic
# identities — no component access required.
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

# ============================================================
# Volume form scalar extraction and VectorValue ↔ 1-form bridge
# ============================================================

"""
    vol_coeff(a::CellField)

Extracts the single scalar coefficient of a top-degree D-form CellField.
Returns an `OperationCellField` (scalar-valued) suitable for use in
`∫(vol_coeff(ω ∧ ⋆η)) * dΩ` with Gridap's standard measure.
"""
vol_coeff(a::CellField) = Operation(vol_coeff)(a)

"""
    to_1form(v::CellField)

Pointwise conversion of a VectorValue{D}-valued CellField (e.g., a Nédélec
FEFunction) to a DifferentialFormValue{1,D}-valued OperationCellField.
The result can be used with the value-level operators (∧, ⋆, ι).
Note: for `exterior_derivative`, build the cell field from `DifferentialForm`
objects with accessible component fields instead.
"""
to_1form(v::CellField) = Operation(to_1form)(v)

"""
    from_1form(ω::CellField)

Pointwise conversion of a DifferentialFormValue{1,D}-valued CellField to a
VectorValue{D}-valued OperationCellField.
"""
from_1form(ω::CellField) = Operation(from_1form)(ω)
