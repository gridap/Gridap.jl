
# Changelog

This file contains the most significant changes to the library as well as upgrade notes for users. For a more detailed list of changes, please refer to the [complete changelog](https://github.com/gridap/Gridap.jl/blob/master/NEWS.md).

## Upgrading to v0.20

This is a major release with breaking changes, mainly in the `Polynomials` and `ReferenceFEs` modules. The main goal of this release is to provide a more flexible and extensible framework for polynomial bases and reference finite elements based on moments, as a stepping stone for the implementation of more complicated elements (e.g. MTW, etc...).

### New features

- **Polynomials**: new `Chebyshev` and `Bernstein` polynomial types; new `CartProdPolyBasis` and `CompWiseTensorPolyBasis` generalising the old monomial bases; new high-level factory `FEEC_poly_basis` covering all spaces in the Periodic Table of the Finite Elements.
- **ReferenceFEs**: new `ReferenceFE(F, r, k[, T]; kwargs...)` constructor using FEEC notation; new `MomentBasedReferenceFE` factory unifying moment-based element construction; new elements `modal_lagrangian`, `modal_serendipity`, `CrouzeixRaviart`, and `nedelec2`; unified Piola map hierarchy (`Pullbacks.jl`) and geometric decomposition API.
- **TensorValues**: new `SkewSymTensorValue` and `HighOrderTensorValue` types; new operations `congruent_prod`, `contracted_product`, and `component_basis`.
- **Default bases and quadratures**: simplex elements (RT, BDM, Nédélec) now use Bernstein bases by default; n-cube elements use Legendre.

### Breaking changes

- **`Nedelec` → `Nedelec{kind}`**: replace `Nedelec()` with `nedelec` (= `Nedelec{1}()`) or `nedelec2` (= `Nedelec{2}()`).
- **`DivConforming` abstract type removed**: `RaviartThomas` and `BDM` now directly subtype `ReferenceFEName`.
- **`CellFE` constructor**: trailing positional args are now keyword args.
- **`SignFlipMap`, `TransformRTDofBasis`, `TransformNedelecDofBasis` removed**: handled internally by the new `Pullbacks.jl` pipeline.
- **`ReferenceFE(name, args...; kwargs...)`** now returns a tuple `(name, args, kwargs)` instead of directly constructing a RefFE.
- **Default polynomial bases changed**: internal representations differ from v0.19 — affects serialised data and exact numerical values.
- **`ThirdOrderTensorValue`** is now a type alias for `HighOrderTensorValue`. Code relying on type identity may break.
- **`MacroReferenceFE` constructor** gains a 3rd positional argument `Name`.
- **`Conformity(reffe, sym::Symbol)` per-element overloads removed**: custom elements must implement `valid_conformity_symbols` instead.
- **`FineToCoarseReferenceFEs.jl` no longer loaded**.
- Old stubs `apply`, `FETerm`, `AffineFETerm`, `LinearFETerm`, `FESource`, `FEEnergy` removed from the `Gridap` namespace.

### Deprecations

| Deprecated | Replacement |
|---|---|
| `MonomialBasis{D}(args...)` | `MonomialBasis(Val(D), args...)` |
| `JacobiPolynomialBasis{D}(args...)` | `LegendreBasis(Val(D), args...)` |
| `JacobiPolynomial` | `Legendre` |
| `QGradMonomialBasis{D}(T, order)` | `FEEC_poly_basis(Val(D), T, order+1, 1, :Q⁻)` |
| `QCurlGradMonomialBasis{D}(T, order)` | `FEEC_poly_basis(Val(D), T, order+1, D-1, :Q⁻; rotate_90=(D==2))` |
| `PCurlGradMonomialBasis{D}(T, order)` | `FEEC_poly_basis(Val(D), T, order+1, D-1, :P⁻; rotate_90=(D==2))` |
| `NedelecPreBasisOnSimplex{D}(order)` | `NedelecPolyBasisOnSimplex(Val(D), Float64, order)` |
| `num_terms(f)` | `length(f)` |
| `return_type(::PolynomialBasis)` | `value_type(::PolynomialBasis)` |
| `n_components` | `num_components` |
| `get_cell_shapefuns(model, cell_reffe, conf)` | `get_cell_shapefuns_and_dof_basis(...)[1]` |
| `get_cell_dof_basis(model, cell_reffe, conf)` | `get_cell_shapefuns_and_dof_basis(...)[2]` |
| `get_free_values(space)` | `get_free_dof_values(space)` |
| `get_dirichlet_values(space)` | `get_dirichlet_dof_values(space)` |
