
"""
    get_cell_metric(model::DiscreteModel)

The metric tensor `g` of the coordinates the cells of `model` map to, either

- `LinearAlgebra.I` the identity metric, or
- a physical cell field `Point{D} -> SymTensorValue{D}`, where `D` is a `Int` or `Cov`.

Currently, `g` needs to be invertible, which requires that `model`'s embedded
dimension equal it's manifold dimension, or that `g` is `LinearAlgebra.I`.

After change of coordinates, the reference volume element of integrals is
`meas(Jt)·√|det g|`: `g` does not take the reference element -> chart element
into account.
"""
get_cell_metric(model::DiscreteModel) = LinearAlgebra.I

"""
    get_cell_inv_metric(model::DiscreteModel)

The inverse `g⁻¹` of the metric tensor, in the shapes of
[`get_cell_metric`](@ref). `g⁻¹` has indices of variance `Con` or `Int`.
"""
get_cell_inv_metric(model::DiscreteModel) = _inv_cell_metric(get_cell_metric(model))

_inv_cell_metric(::UniformScaling) = LinearAlgebra.I
_inv_cell_metric(cell_metric::AbstractArray) = lazy_map(Operation(inv), cell_metric)
