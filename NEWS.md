# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.18.6] - 2024-08-29

### Fixed

- Improved performance of PR[#967](https://github.com/gridap/Gridap.jl/pull/967). Along the way, opened the door to Triangulations of different type in SkeletonTriangulation. Since PR[#1026](https://github.com/gridap/Gridap.jl/pull/1026).

## [0.18.5] - 2024-08-28 

### Changed

- Misc changes required to support facet integration on non-conforming meshes. These changes do not involve methods of the public API. Since PR[#967](https://github.com/gridap/Gridap.jl/pull/967)

## [0.18.4] - 2024-08-09

### Changed

- Added WriteVTK kwargs to control the output encoding for vtk files. Since PR[#1016](https://github.com/gridap/Gridap.jl/pull/1016).

### Fixed

- Passing `kwargs` from `refine` to `simplexify` functions in Adaptivity. Since PR[#1015](https://github.com/gridap/Gridap.jl/pull/1015).
- Fixed `interpolate` for `ZeroMeanFESpace`. Since PR[#1020](https://github.com/gridap/Gridap.jl/pull/1020).
- Fixed `gather_free_and_dirichlet_values!` for `FESpaceWithConstantFixed`. Since PR[#1020](https://github.com/gridap/Gridap.jl/pull/1020).

## [0.18.3] - 2024-07-11

### Added

- Added more features to Adaptivity. Notably: 3D uniform edge-based refinement for HEX and TETs. Barycentric refinement for simplices. Simplexify as a new refinement strategy. Since PR[#1013](https://github.com/gridap/Gridap.jl/pull/1013).

- Define `GeneralPolytope` that represents general polytopes in 2 and 3 dimensions. Since PR[#1006](https://github.com/gridap/Gridap.jl/pull/1006).

### Changed

- Added local preferences for Gridap through `Preferences.jl`. The macro `@check` can now be activated/deactivated by using the local preference `execution_mode`. Since PR[#1014](https://github.com/gridap/Gridap.jl/pull/1014).

## [0.18.2] - 2024-05-02

### Fixed

- Bugfix in `get_face_dofs` for Nedelec GenericRefFE. Since PR[#1005](https://github.com/gridap/Gridap.jl/pull/1005).
- Ensure deterministic behavior for matrix assembly involving multiple domains. Since PR[#1004](https://github.com/gridap/Gridap.jl/pull/1004).

## [0.18.1] - 2024-04-12

### Changed

- Changed the sign of the residual in `TransientLinearFEOperator` to align with the conventions of `AffineFEOperator`. Since PR[#996](https://github.com/gridap/Gridap.jl/pull/996).

### Fixed

- Bugfix in `restrict_to_field` for `BlockMultiFieldStyle`. Since PR[#993](https://github.com/gridap/Gridap.jl/pull/993).

## [0.18.0] - 2024-04-12

### Breaking

- ODE module extensive refactor. Breaking changes! See docs and PR for details. Since PR[965](https://github.com/gridap/Gridap.jl/pull/965).
- Fixed name clash with `Statistics.mean`. Since PR[#988](https://github.com/gridap/Gridap.jl/pull/988).
- Deprecated `SubVector` in favor of Julia's `view`. Since PR[#989](https://github.com/gridap/Gridap.jl/pull/989).

### Added

- Added some missing API methods to `Assemblers` and `MultiField`. Since PR[#985](https://github.com/gridap/Gridap.jl/pull/985).

### Fixed

- Fix when evaluating `\circ` operator with `CellState`. Since PR[#987](https://github.com/gridap/Gridap.jl/pull/987).

## [0.17.23] - 2024-01-28

### Changed

- Changed how `symbolic_loop_matrix_vector!` loop works. Now it also takes account vector entries touched from matvecdata. Since PR[#975](https://github.com/gridap/Gridap.jl/pull/975).

## [0.17.22] - 2024-01-12

### Added

- Added `get_dof_to_node` and `get_dof_to_comp` for `LagrangianDofBasis`. Since PR[#964](https://github.com/gridap/Gridap.jl/pull/964).

### Changed

- Changed how `allocate_vector` works. Now it only allocates, instead of allocating+initialising to zero. Since PR[#963](https://github.com/gridap/Gridap.jl/pull/963).

### Fixed

- Fixed issue where `FineToCoarseField` would not work for vector-valued fields. Since PR[#970](https://github.com/gridap/Gridap.jl/pull/970).
- Fixed issue where `BlockSparseMatrixAssembler` would not work correctly when permuting variables. Since PR[#971](https://github.com/gridap/Gridap.jl/pull/971).

## [0.17.21] - 2023-12-04

### Added

- Implemented real/imag for VectorValues
- Explicit Runge-Kutta ODE Solvers. Since PR [#952](https://github.com/gridap/Gridap.jl/pull/952)
- Improved the methods `allocate_in_range` and `allocate_in_domain` with support for `BlockArrays` and distributed arrays. Since PR[#960](https://github.com/gridap/Gridap.jl/pull/960).

### Fixed

- `BlockMultiFieldStyle` available for `TransientMultiFieldFESpaces` since PR [#946](https://github.com/gridap/Gridap.jl/pull/946).
- When creating `DiscreteModelPortions`, some of the `FaceLabeling` arrays were being aliased. This caused issues when adding tags to distributed models in debug mode. Since PR [#956](https://github.com/gridap/Gridap.jl/pull/956).
- Function `add_entry!` was inconsistent for `AbstractMatrix` and `AbstractSparseMatrix`. Since PR[#959](https://github.com/gridap/Gridap.jl/pull/959).

## [0.17.20] - 2023-10-01

### Added

- Block assembly now generalised to work with `AbstractBlockArrays`, to include changes in GridapDistributed. Since PR [939](https://github.com/gridap/Gridap.jl/pull/939).
- Implici-Explicit Runge-Kutta ODE solvers. Since PR [#919](https://github.com/gridap/Gridap.jl/pull/919).

### Fixed

- Using Broadcasting(\circ) instead of \circ in one of the lazy_maps used to transform a coarse field into a fine field. Since PR [#938](https://github.com/gridap/Gridap.jl/pull/938).
- Better infinite norm computation in `Algebra._check_convergence`. Now works for any `AbstractArray` type, including `PVector`. Since PR [#940](https://github.com/gridap/Gridap.jl/pull/940).
- Updated Runge-Kutta solver. Since PR [#919](https://github.com/gridap/Gridap.jl/pull/919).

## [0.17.19] - 2023-08-23

### Fixed

- Reimplemented `DomainStyle` for `CellQuadrature` to fix breaking low-level Poisson tutorial. Since PR [#937](https://github.com/gridap/Gridap.jl/pull/937).

## [0.17.18] - 2023-08-15

### Added

- Jacobi polynomial bases. Since PR [#896](https://github.com/gridap/Gridap.jl/pull/896).
- Replaced newest vertex bisection mesh adaptation in
  `src/Geometry/NewestVertexBisection.jl` with appropriate changes to
  `src/Adaptivity/EdgeBasedRefinement.jl`. Since PR
  [#901](https://github.com/gridap/Gridap.jl/pull/901).
- When refining `DiscreteModels`, the `FaceLabeling` of the resulting `AdaptedDiscreteModel` will now correctly inhering the tags of the parent model. This has been made possible by the addition of the method `get_d_to_face_to_parent_face`. Since PR[#886](https://github.com/gridap/Gridap.jl/pull/886).
- Added support for mixed adaptivity (i.e coarsening and refining), as well as non-conforming adaptivity. Since PR[#886](https://github.com/gridap/Gridap.jl/pull/886).
- Added support for block assembly of FE systems. Two new types `BlockMultiFieldStyle` and `BlockSparseMatrixAssemblers` have been added. Since PR[#915](https://github.com/gridap/Gridap.jl/pull/915).

### Changed

- The API of `CellQuadrature` has now both data and integration domain styles as keyword arguments. Old signatures are deprecated. Since PR [#885](https://github.com/gridap/Gridap.jl/pull/885).

### Fixed

- ODE operators cache linear system at initial time or the time stored by the operator. Before, the linear system was cached at time `t = 0.0`, which cannot be done if the operator is not well-defined at `t = 0.0`. Since PR [#891](https://github.com/gridap/Gridap.jl/pull/891).
- Fixed the method `get_normal_vector` for `AdaptedTriangulation`. The method `get_facet_normal`
  was using default, it's now using the spetialized implementation for the underlying triangulation type.
  Since PR [#884](https://github.com/gridap/Gridap.jl/pull/884).
- Fixed `cell_dof_ids` for the case of vectorial `ConstantFESpace`. Since PR [#888](https://github.com/gridap/Gridap.jl/pull/888)
- Fixed generation of Modal C0 bases for Julia 1.9. Since PR [#918](https://github.com/gridap/Gridap.jl/pull/918).
- Fixed some edge cases for `change_domain` between `AdaptedTriangulations` where inneficient coordinate transformations would be applied between physical and reference domains. Since PR[#886](https://github.com/gridap/Gridap.jl/pull/886).
- Fixed: Domain limits can now be of any type (notably, floats) when refining `CartesianDiscreteModels`. Since PR[#886](https://github.com/gridap/Gridap.jl/pull/886).

## [0.17.17] - 2023-02-24

### Added

- Implemented `RefinementRule` and `AdaptivityGlue`, which encode the mapping between the old and new cells in an adapted model. Since PR [#838](https://github.com/gridap/Gridap.jl/pull/838).
- Implemented `AdaptedDiscreteModel` and `AdaptedTriangulation`, representing respectively a model and triangulation produced by adapting a parent model. This types mostly wrap around `DiscreteModel` and `Triangulation`, with some added features necessary to keep track of the adaptive hierarchy. Since PR [#838](https://github.com/gridap/Gridap.jl/pull/838).
- Implemented tools to be able to transfer `CellDatum`s back and forth between parent and child grids. These include changes to `change_domain` (which now takes the source `Triangulation` as argument) and a new type of `Measure` called `CompositeMeasure`, which allows the integration `Triangulation` to be different from the `Triangulation` of the resulting `DomainContribution`. To accommodate `CompositeMeasure`, `Measure` has been made abstract type and a `GenericMeasure` has been created to replace the old type. Since PR [#838](https://github.com/gridap/Gridap.jl/pull/838).
- For the fine-to-coarse transfer of `CellField`s, the new `FineToCoarseField` has been implemented. This new structure bundles several fields defined on the fine mesh to create a single field on the coarse mesh. To enable fast interpolation of this type of field, we have also implemented `FineToCoarseReferenceFE` and `FineToCoarseDofBasis`. Since PR [#838](https://github.com/gridap/Gridap.jl/pull/838).
- Implemented `CompositeQuadrature`, a quadrature for a cell that has been refined using a `RefinementRule`. Since PR [#838](https://github.com/gridap/Gridap.jl/pull/838).
- Implemented simple refinement strategies for Cartesian discrete models in 2&3D as well as Unstructured discrete models in 2D. The latter is implemented by red-green refinement. Since PR [#838](https://github.com/gridap/Gridap.jl/pull/838).
- Added optimization when calling `unique` for a `CompressedArray`. Since PR [#838](https://github.com/gridap/Gridap.jl/pull/838).
- Added support for changing domain between adapted triangulations in cases where the target triangulation is a view, a `BoundaryTriangulation` or a `SkeletonTriangulation`. Since PR [#868](https://github.com/gridap/Gridap.jl/pull/868).

### Fixed

- Using broadcasting through in `ODESolver` vector operations. Since PR [#860](https://github.com/gridap/Gridap.jl/pull/860)
- Fixes to `array_cache(a::Table)`: Now does not use the `zero(T,N)` function, but instead creates new empty vector using the general allocator `Vector{T}(undef,N)`. This allows `Table` to work with complex composite types which don't have an easy `zero(T)` function defined. Since PR [#838](https://github.com/gridap/Gridap.jl/pull/838).
- Added `get_metadata` to all the instances of `ReferenceFE`. This makes the abstract type more consistent, which is necessary for the new type `FineToCoarseReferenceFE`. Since PR [#838](https://github.com/gridap/Gridap.jl/pull/838).
- `ConstantFESpace` is now properly exported. Since PR [#872](https://github.com/gridap/Gridap.jl/pull/872).

## [0.17.16] - 2022-12-22

### Fixed

- Fixed a bug in function `collect_cell_matrix_and_vector`. Since PR [#849](https://github.com/gridap/Gridap.jl/pull/849)

## [0.17.15] - 2022-11-10

### Added

- `lastindex` for `MultiValue`s for consistent usage of `[end]` as per `length`. Since PR [#834](https://github.com/gridap/Gridap.jl/pull/834)
- BDM (Brezzi-Douglas-Marini) ReferenceFEs in PR [#823](https://github.com/gridap/Gridap.jl/pull/823)
- Implemented `ConstantFESpace`. This space allows, e.g., to impose the mean value of the pressure via an additional unknown/equation (a Lagrange multiplier). Since PR [#836](https://github.com/gridap/Gridap.jl/pull/836)
- Methods to construct `DiracDelta` at generic `Point`(s) in the domain. Since PR [#837](https://github.com/gridap/Gridap.jl/pull/837)
- some key missing `lastindex` and `end` tests belonging to PR [#834]. Since PR [#837](https://github.com/gridap/Gridap.jl/pull/837)
- `AbstractVector` support to weights and `Point`s of GenericQuadrature, in particular for `FillArray` usage. Since PR [#839](https://github.com/gridap/Gridap.jl/pull/839)
- Using the above generalization of `GenericQuadrature`, made the weights of quadrature of `DiracDelta` for generic points to be a `Fill` vector. Since PR [#839](https://github.com/gridap/Gridap.jl/pull/839)
- Optional keyword argument `T` to select the floating point precision for the numerical quadrature pipeline involving `Measure`, `CellQuadrature` and `Quadrature` functions. `T` defaults to `Float64` when not specified, so it doesn't involve any breaking changes in public API. Since PR [#840](https://github.com/gridap/Gridap.jl/pull/840)

## [0.17.14] - 2022-07-29

### Added

- Functionality to take gradient of functional involving integration (`DomainContribution`) over Skeleton faces, with respect to the degrees-of-freedom of `FEFunction`. The interface remains the same - `gradient(f,uh)`. Since PR [#797](https://github.com/gridap/Gridap.jl/pull/797)
- Extended the `MultiField` functional gradient (with respect to degrees-of-freedom of `MultiFieldFEFunction`) to functionals involving Skeleton integration. The interface remains the same `gradient(f,xh)`. Since PR [#799](https://github.com/gridap/Gridap.jl/pull/799)
- Functionality to take jacobian of functional involving integration (`DomainContribution`) over Skeleton faces (obtained from testing bilinear form with the whole set of test `fe_basis`), with respect to the degrees-of-freedom of `FEFunction`. The interface remains the same - `jacobian(f,uh)`. Since PR [#803](https://github.com/gridap/Gridap.jl/pull/803)
- A dummy tag for ForwardDiff configs being constructed in Gridap at `src/Arrays/Autodiff.jl` to fix issue [#805](https://github.com/gridap/Gridap.jl/issues/805). Since PR [#806](https://github.com/gridap/Gridap.jl/pull/806)

### Fixed

- The behavior of `gradient` for functionals involving operations of `CellFields` inside `mean` and `jump` terms of Skeleton Integration terms. Since PR [#800](https://github.com/gridap/Gridap.jl/pull/800)
- The behavior of `SkeletonCellFieldPair` at the Boundary integration terms. Since PR [#800](https://github.com/gridap/Gridap.jl/pull/800)
- `push_normal` for rectangular Jacobians. Since PR[#809](https://github.com/gridap/Gridap.jl/pull/809)
- Nedelec FEs with triangles. Since PR [#813](https://github.com/gridap/Gridap.jl/pull/813)
- `Gridap.jacobian` for Skeleton integration terms, in the case of general test basis where dim not same as trial basis. Since PR [#815](https://github.com/gridap/Gridap.jl/pull/815)

## [0.17.13] - 2022-05-31

### Added

- `KeyToValMap` lazy map that dynamically creates a `Dict` with the outputs of a function over an array of inputs. Since PR [#801](https://github.com/gridap/Gridap.jl/pull/801)
- `MappedDiscreteModel` and `MappedGrid`, which are geometrical models with one extra geo map in the physical space. Since PR [#801](https://github.com/gridap/Gridap.jl/pull/801)
- `GridWithFEMap`, which has a geometrical map defined by a FE function. Since PR [#801](https://github.com/gridap/Gridap.jl/pull/801)
- Vertex bisection algorithm for refinement of triangular meshes in 2D. Since PR [#733](https://github.com/gridap/Gridap.jl/pull/733)
- Generalized-α method for 1st order ODEs. Since PR [#781](https://github.com/gridap/Gridap.jl/pull/781)
- Implemented (generalised) ModalC0 Polynomial bases and reference FEs. Since PR [#777](https://github.com/gridap/Gridap.jl/pull/777)
- Serendipity reference FEs for any dimension and order. Since PR [#777](https://github.com/gridap/Gridap.jl/pull/777)

## [0.17.12] - 2022-03-24

### Fixed

- Bugs related with 2nd order derivatives. Since PR [#769](https://github.com/gridap/Gridap.jl/pull/769) and [#771](https://github.com/gridap/Gridap.jl/pull/771).

## [0.17.11] - 2022-03-16

### Fixed

- Miscellaneous bugs associated with empty meshes. Since PR [#765](https://github.com/gridap/Gridap.jl/pull/765) and [#766](https://github.com/gridap/Gridap.jl/pull/766).

## [0.17.10] - 2022-03-10

### Added

- Implemented addition/subtraction between a `TensorValue` and a `SymTensorValue`. Since PR [#755](https://github.com/gridap/Gridap.jl/pull/755).
- Nédélec elements on simplices (lowest order interpolation for the moment). Since PR [#758](https://github.com/gridap/Gridap.jl/pull/758).

### Fixed

- Restrict to the active model the evaluation of a FE function at arbitrary points. Since PR [#752](https://github.com/gridap/Gridap.jl/pull/752).
- Renamed `GridapODEs` to `ODEs`. Since PR [#756](https://github.com/gridap/Gridap.jl/pull/756)

## [0.17.9] - 2022-02-28

### Added

- Extra support for dual number propagation. Since PR [#653](https://github.com/gridap/Gridap.jl/pull/653)
- Added ODE tools to Gridap. [GridapODEs](https://github.com/gridap/GridapODEs.jl) has been merged into Gridap. Sice PR[#748](https://github.com/gridap/Gridap.jl/pull/748)

## [0.17.8] - 2022-02-14

### Added

- Support for periodic conditions in `CartesianDiscreteModel`s built from `cmin`, `cmax`. Since PR [#738](https://github.com/gridap/Gridap.jl/pull/738).
- Function `savepvd`. Since PR [#736](https://github.com/gridap/Gridap.jl/pull/736).
- Support for passing `nothing` as the first argument of the `CartesianDiscreteModel` constructor to mimic the API of `GridapDistributed`. Since PR [#735](https://github.com/gridap/Gridap.jl/pull/735).

### Fixed

- Bug in `SymTensorValue`. Since PR [#740](https://github.com/gridap/Gridap.jl/pull/740).
- Bug in interpolation at arbitrary points. Since PR [#729](https://github.com/gridap/Gridap.jl/pull/729).
- Bug in `reference_grid(TET)`. Since PR [#734](https://github.com/gridap/Gridap.jl/pull/734).
- Bug in computation of Hessian with AD for multi-field. Since PR [#731](https://github.com/gridap/Gridap.jl/pull/731).

## [0.17.7]- 2021-12-03

### Fixed

- Improving compile times by removing unnecessary `@inline` statements. Since PR [#726](https://github.com/gridap/Gridap.jl/pull/726).

### Added

- CellFE constructor now gets optional arguments and pass them down. Since PR [#728](https://github.com/gridap/Gridap.jl/pull/728).

## [0.17.6] - 2021-12-01

### Added

- Implemented `Base.unaliascopy(A::SubVector)`. Since PR [#715](https://github.com/gridap/Gridap.jl/pull/715).

### Fixed

- Bugfix in `Base.view(glue::FaceToFaceGlue,ids::AbstractArray)`. Since PR [#724](https://github.com/gridap/Gridap.jl/pull/724).

## [0.17.5] - 2021-11-08

### Added

- Hiding the creation of `FESolver` and `LinearFESolver` from user code. Since PR [#705](https://github.com/gridap/Gridap.jl/pull/705).

## [0.17.4] - 2021-10-27

### Fixed

- Using implementation of `pvtk_grid` provided in WriteVTK. Since PR [#699](https://github.com/gridap/Gridap.jl/pull/699).

## [0.17.3] - 2021-10-27

### Fixed

- Adding a newer version of WriteVTK in the [compat] section. Since PR [#698](https://github.com/gridap/Gridap.jl/pull/698).

## [0.17.2] - 2021-10-26

### Fixed

- AD with multi-field residuals with different num dofs per field. Since PR [#687](https://github.com/gridap/Gridap.jl/pull/687).

## [0.17.1] - 2021-10-26

### Fixed

- Laplacian `Δ` operator on unstructured linear grids for quantities defined in the reference space (i.e. shape functions in standard FEM). Since PR [#691](https://github.com/gridap/Gridap.jl/pull/691).
- Laplacian `Δ` operator on triangulations using `GridView` (e.g., when interpolating functions in a sub-domain or on the boundary). Since PR [#691](https://github.com/gridap/Gridap.jl/pull/691).
- Fixed typo in , function `solve! of `LinearSolvers.jl`. Since PR [#692](https://github.com/gridap/Gridap.jl/pull/692).

## [0.17.0] - 2021-10-22

### Added

- Aliases `Interior`, `Boundary`, `Skeleton`, and `Interface` for the `Triangulation`, `BoundaryTriangulation`, `SkeletonTriangulation`, and `InterfaceTriangulation` constructors. Since PR [#662](https://github.com/gridap/Gridap.jl/pull/662).
- Function `create_pvtk_file` for exporting results in `pvtu` format. Since PR [#685](https://github.com/gridap/Gridap.jl/pull/685).

### Changed

- Major refactoring in the `Triangulation` interface to properly support the solution of PDEs defined on domains of different dimension. The major change from the user perspective is that `Triangulation` objects can be used both to integrate the weak form (as before) but also to define FE spaces (except for unfitted triangulations obviously). It is still possible to define FE spaces from `DiscreteModels`, but it is safer and more idiomatic (closer to the math notation) to use `Triangulation` objects from now on. Since PR [#662](https://github.com/gridap/Gridap.jl/pull/662).
- Changes in assembly interface to allow optimization when assembling matrices and vectors simultaneously. Since PR [#685](https://github.com/gridap/Gridap.jl/pull/685).

### Removed

- `BoundaryDiscreteModel`, `RestrictedDiscreteMdeol`, `RestrictedTriangulation`, `TriangulationStyle`, `BackgroundTriangulation`, `SubTriangulation`, `get_cell_to_bgcell`, `get_cell_ref_map`, `get_background_triangulation`, and `have_compatible_domains`. Since PR [#662](https://github.com/gridap/Gridap.jl/pull/662).
- Functions `scale_entries!` and `fill_entries!`. Replaced by Julia functions `LinearAlgebra.rmul!` and `LinearAlgebra.fillstored!`. Since PR [#680](https://github.com/gridap/Gridap.jl/pull/680).

## [0.16.5] - 2021-09-08

### Added

- Implemented DIV operator for FE functions in RT space. Since PR [#650](https://github.com/gridap/Gridap.jl/pull/650).
- `GenericAssemblyStrategy`. Since PR [#655](https://github.com/gridap/Gridap.jl/pull/655).
- Additional high level API assembly functions. Since PR [#652](https://github.com/gridap/Gridap.jl/pull/652).

### Fixed

- Bug related with the release of ChainRulesCore version 1.3.1. Since [#654](https://github.com/gridap/Gridap.jl/pull/654).
- Inheritance relationship for DiscreteModelPortion. Since PR [#645](https://github.com/gridap/Gridap.jl/pull/645).
- Optimization to RT FEs. Since PR [#638](https://github.com/gridap/Gridap.jl/pull/638).
- Bug in boundary discrete model. Since PR [#651](https://github.com/gridap/Gridap.jl/pull/651).

## [0.16.4] - 2021-08-17

### Added

- Trait to CellQuadrature to support the evaluation of integrals on the reference domain. Since PR [#636](https://github.com/gridap/Gridap.jl/pull/636).
- Type `Interpolable` allowing to interpolate data from two independent meshes. Since PR [#632](https://github.com/gridap/Gridap.jl/pull/632).

## [0.16.3] - 2021-06-28

### Fixed

- Deactivating optimizations related with `MemoArray` since they are not reliable. Since PR [#624](https://github.com/gridap/Gridap.jl/pull/624).
- Bug related with `ArrayBlock`. Since PR [#623](https://github.com/gridap/Gridap.jl/pull/623).

## [0.16.2] - 2021-06-21

### Fixed

- Bug related with boundary integration caused by some optimization introduced in v0.16. Fixed via PR [#616](https://github.com/gridap/Gridap.jl/pull/616).

## [0.16.1] - 2021-06-04

### Fixed

- Bug for 1st order FE spaces in combination of 1st order models with periodic BCs. Since PR [#611](https://github.com/gridap/Gridap.jl/pull/611).

## [0.16.0] - 2021-06-04

### Added

- User API to select specific quadrature rules. Since PR [#578](https://github.com/gridap/Gridap.jl/pull/578).
- Experimental support for mixed dimensional PDEs. Since PR [#567](https://github.com/gridap/Gridap.jl/pull/567).
- Added `get_cell_dof_basis(model,cell_reffes,::Conformity)` and `get_cell_shapefuns(model,cell_reffes,::Conformity)`. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).
- Implemented `get_cell_dof_basis` and `get_cell_shapefuns` for global RT FE spaces in a new file `DivConformingFESpaces.jl`. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).
- Added support to allow evaluation of FE functions at arbitrary points. Since PR [#523](https://github.com/gridap/Gridap.jl/pull/523).
- Implemented `compute_cell_points_from_vector_of_points` to build `CellPoint` from a vector of points. Since PR [#523](https://github.com/gridap/Gridap.jl/pull/523).

### Changed

- Major refactoring in the handling of blocks (e.g. in multi-field and skeleton terms). The new code follows a much more simple approach based in the new type `ArrayBlock`. Since PR [#583](https://github.com/gridap/Gridap.jl/pull/583).
- The default quadrature rule for tets has changed. Since PR [#578](https://github.com/gridap/Gridap.jl/pull/578).
- Refactoring in `SparseMatrixAssembler` to make it more extensible and efficient. Since PR [#568](https://github.com/gridap/Gridap.jl/pull/568).
- Renamed `get_free_values` -> `get_free_dof_values`. Since PR [#567](https://github.com/gridap/Gridap.jl/pull/567).
- Renamed `get_dirichlet_values` -> `get_dirichlet_dof_values`. Since PR [#606](https://github.com/gridap/Gridap.jl/pull/606).
- Renamed `object` -> `value` the variable in `ConstantField`. Since PR [#606](https://github.com/gridap/Gridap.jl/pull/606).
- Miscellaneous changes in the FE assembly to allow the solution of mixed dimensional problems. Since PR [#567](https://github.com/gridap/Gridap.jl/pull/567).
- Renamed `get_cell_shapefuns` by `get_fe_basis`. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).
- Renamed `get_cell_shapefuns_trial` by `get_trial_fe_basis`. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).
- Renamed `get_cell_dof_basis` by `get_fe_dof_basis`. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).
- Removed `conformity` optional keyword argument from `FESpace(::DiscreteModel,::CellFE; kwargs...)` constructor. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).
- Replaced `CellFE(::AbstractArray{<:Field},::AbstractArray{<:ReferenceFE})` by `CellFE(::DiscreteModel,::AbstractArray{<:ReferenceFE},::Conformity)`. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).

### Removed

- All code associated with with `BlockArrayCoo`. Since PR [#583](https://github.com/gridap/Gridap.jl/pull/583).
- Module `Gridap.Integration` has been deleted and its contents have been merged into `Gridap.ReferenceFEs` module.
- Types `SparseMatrixCSR` and `SymSparseMatrixCSR` have been moved to the registered package [`SparseMatricesCSR`](https://github.com/gridap/SparseMatricesCSR.jl). To use them simply add `SparseMatricesCSR` into your environment and type `using SparseMatricesCSR`. Since Since PR [#568](https://github.com/gridap/Gridap.jl/pull/568).
- Removed `PushForwardMap` and all code depending upon it. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).

## [0.15.5] - 2021-05-25

### Added

- Differential operators `(∇+k)(u)`, `(∇+k)⋅u`, `(∇+k)×u`, `(∇+k)⊗u`, and `u⊗(∇+k)` for some `u::CellField` and `k::VectorValue`. Since PR [#597](https://github.com/gridap/Gridap.jl/pull/597).
- Definition of `u.*v` between instances of vector-valued `CellField` objects `u` and `v`. Also differential operators `∇.*u` and `(∇+k).*u`. Since PR [#597](https://github.com/gridap/Gridap.jl/pull/597).

## [0.15.4] - 2021-03-29

### Fixed

- Bug in `CartesianDiscreteModel` with periodic boundary conditions that shows up in Julia 1.6 but not in Julia 1.5. Since commit [da005cf](https://github.com/gridap/Gridap.jl/commit/da005cf4cde68617f92d76744e307798ef7e8340).

## [0.15.3] - 2021-03-16

### Added

- `get_cell_map` now returns array of `AffineMap` for linear grids of simplices. Needed to compute Laplacian operator, inverse maps etc. Since PR [#553](https://github.com/gridap/Gridap.jl/pull/553).

### Fixed

- Bug in `print_op_tree`. Since PR [#563](https://github.com/gridap/Gridap.jl/pull/563)

## [0.15.2] - 2021-03-08

### Added

- Method `inverse_map` for `AffineMap`. Since PR [#552](https://github.com/gridap/Gridap.jl/pull/552).
- Method `get_cell_points` for `CellDof`. Since PR [#551](https://github.com/gridap/Gridap.jl/pull/551).
- Evaluation of `MonomialBasis` objects at a single point. Since PR [#550](https://github.com/gridap/Gridap.jl/pull/550).
- `rand` function for `MultiValue` objects. Since PR [#530](https://github.com/gridap/Gridap.jl/pull/530).

### Fixed

- Bug in `return_value` for `Broadcasting(∇∇)`. Since PR [#554](https://github.com/gridap/Gridap.jl/pull/554).
- Bug in `dot` for third order tensors. Since PR [#544](https://github.com/gridap/Gridap.jl/pull/544).

## [0.15.1] - 2021-01-22

### Added

- Added support for Hessian and Laplacian operators. Only implemented for Finite Elements with an `AffineMap`. Since PR [#514](https://github.com/gridap/Gridap.jl/pull/514).

### Fixed

- Bug in `RestrictedDiscreteModel` for periodic boundary conditions. Since PR [#517](https://github.com/gridap/Gridap.jl/pull/517).
- Bug in `sum(a::LazyArray)` when `eltype(a) <: AbstractArray`. Since PR [#513](https://github.com/gridap/Gridap.jl/pull/513).

## [0.15.0] - 2020-12-14

This version is a major (backwards-incompatible) refactoring of the project which is not summarized here for the sake of brevity. Most of the functionality of v0.14.0 is available in v0.15.0, but possibly with a significantly different API. See [here](https://github.com/gridap/Tutorials/compare/v0.14.0...v0.15.0) the changes in the sources of the Gridap Tutorials between versions 0.14.0 and 0.15.0 to effectively see the major changes in the API.

## [0.14.2] - 2020-11-24

### Added

- Added additional tensor operations and new double contraction notation `⋅²`. Implemented a `zero` constructor for `ThirdOrderTensorValues` to allow integration of 3-tensors. Since PR [#415](https://github.com/gridap/Gridap.jl/pull/415/).

### Fixed

- Bug-fix for 32-bit Julia: Replace all occurrences of Int64 by Int. Since PR [#445](https://github.com/gridap/Gridap.jl/pull/445).
- Bug-fix for 32-bit Julia. Using inttype=Int keyword argument for JSON parsing. Since PR [#456](https://github.com/gridap/Gridap.jl/pull/456).

## [0.14.1] - 2020-09-17

### Added

- Added VectorWithEntryInserted and VectorWithEntryRemoved. Since PR [#401](https://github.com/gridap/Gridap.jl/pull/401/).
- Added missing get_constant_approach() getter to FESpaceWithConstantFixed. Since PR [#409](https://github.com/gridap/Gridap.jl/pull/409).

### Deprecated

- The name FESpaceWithLastDofRemoved has been deprecated in favor of its generalization FESpaceWithConstantFixed. Since PR [#396](https://github.com/gridap/Gridap.jl/pull/396) and PR [#404](https://github.com/gridap/Gridap.jl/pull/404).

## [0.14.0] - 2020-08-27

### Removed

- Support for Julia v1.0. Now, the minimum supported is Julia v1.3. Since PR [#376](https://github.com/gridap/Gridap.jl/pull/376/).

### Changed

- Major refactoring associated with the handling of elemental matrices and vectors in multi-field computations and also on the skeleton. Since PR [#376](https://github.com/gridap/Gridap.jl/pull/376/).
- First and second argument switch in `update_state_variables!` in order to have function-first style. Since PR [#376](https://github.com/gridap/Gridap.jl/pull/376/).
- Table struct has been generalized such that data and ptrs arrays can be of an arbitrary type extending AbstractArray. Since PR [#310](https://github.com/gridap/Gridap.jl/pull/310/)
- `interpolate, interpolate!, interpolate_dirichlet...` switched argument order to function first style. For instance `interpolate(u, V)` instead of `interpolate(V, u)`

### Added

- Allowing the construction of an `HomogeneousTrialFESpace` from a `TrialFESpace`. Since PR [#384](https://github.com/gridap/Gridap.jl/pull/384).
- Support for automatic differentiation of residuals and Jacobians in multi-field computations since PR [#383](https://github.com/gridap/Gridap.jl/pull/383/).
- New `FilterKernel` since PR [#379](https://github.com/gridap/Gridap.jl/pull/379/).

### Fixed

- Bug associated with boundary triangulation in 1D discrete models. Since PR [#393](https://github.com/gridap/Gridap.jl/pull/393).

## [0.13.4] - 2020-08-23

### Added

- New `FilteredCellArray` since PR [#372](https://github.com/gridap/Gridap.jl/pull/372/).

## [0.13.3] - 2020-08-12

### Added

- `Visualization.visualization_data` function that makes it easier to bring fields into
  visualization library friendly formats. Since PR [#354](https://github.com/gridap/Gridap.jl/pull/354).
- Gradient of a product binary operation (`*`) between a scalar and a field. Since PR [#340](https://github.com/gridap/Gridap.jl/pull/340).

## [0.13.2] - 2020-07-31

### Added

- Automatic differentiation of the Jacobian from a given residual and the Jacobian and the residual from a given energy. Not working at this moment on the Skeleton nor for multi-field (WIP), but yes for other cases.
  Now, the user can omit `jac` from `FETerm(res,jac,trian,quad)`, i.e. `FETerm(res,trian,quad)` and the Jacobian will be automatically generated. In addition, the user can write `FEEnergy(ener,trian,quad)` for a given `ener(uh)` function
  and the residual and the Jacobian will be automatically generated. Since PR [#338](https://github.com/gridap/Gridap.jl/pull/338/).

## [0.13.1] - 2020-07-24

### Fixed

- Bugs associated with the degenerated case of 0-length arrays. Since PR [#331](https://github.com/gridap/Gridap.jl/pull/331/) and [#332](https://github.com/gridap/Gridap.jl/pull/332/).

## [0.13.0] - 2020-07-23

### Added

- Automatic differentiation for symmetric gradient, i.e. `ε(u)` for a given vector-valued function `u`. Since PR [#327](https://github.com/gridap/Gridap.jl/pull/327/).
- Added missing SparseMatrixAssembler constructor for MultiFieldFESpaces. Since PR [#320](https://github.com/gridap/Gridap.jl/pull/320/).
- kw-argument `space` to `LagrangianRefFE` constructor in order to select the type of underlying polynomial space, i.e., `:Q`, `:S`, or `:P`. Since PR [#321](https://github.com/gridap/Gridap.jl/pull/321).

### Changed

- The meaning of `inward/outward` has slightly changed for `SkeletonCellBasis` objects. Now, by accessing to these properties a `ReducedSkeletonCellBasis` is returned, which allows to use the result in a more flexible way (in particular, the result can be used in a similar way than the result of `jump` or `mean`). Since PR [#317](https://github.com/gridap/Gridap.jl/pull/317).
- Major refactoring in `ReferenceFEs` module. Since PR [#319](https://github.com/gridap/Gridap.jl/pull/319) and [#321](https://github.com/gridap/Gridap.jl/pull/321). In particular:
  - `NodalReferenceFE` has been replaced by a new abstract type `LagrangianRefFE`.
  - `GenericNodalCartesianRefFE` has been replaced by `GenericLagrangianRefFE`.

### Removed

- Removals associated with the `ReferenceFEs` refactoring in PR [#319](https://github.com/gridap/Gridap.jl/pull/319):
  - Removed `QDiscRefFE` constructor. Use a standard `LagrangianRefFE` and `L2Conformity` instead.
  - Removed `PDiscRefFE` constructor. Use `LagrangianRefFE` constructor with the kw-argument `space=:P`.
  - Removed `CDLagrangianRefFE` constructor. Use a standard `LagrangianRefFE` and `CDConformity` instead.
  - Removed fields `face_own_dofs` and `face_own_dof_permutations` from `GenericRefFE`.
  - Removed struct `DiscRefFE`.

### Fixed

- Better handling of FE terms defined on empty triangulations. Since PR [#329](https://github.com/gridap/Gridap.jl/pull/329).
- Replaced `+=` by `add_entry!`. Since PR [#316](https://github.com/gridap/Gridap.jl/pull/316).
- Minor fix to let Vtk.jl support changes in Vtk 1.7.X versus 1.6.X. Since PR [#324](https://github.com/gridap/Gridap.jl/pull/324).

## [0.12.0] - 2020-07-07

### Added

- Added `SkeletonTriangulation` constructor in order to integrate, where a given interpolation is discontinuous. Since PR [#304](https://github.com/gridap/Gridap.jl/pull/304).
- New `ConformingFESpace` constructor. Since PR [#293](https://github.com/gridap/Gridap.jl/pull/293).
- Added `QDiscRefFE` constructor for `DiscRefFE`. Since PR [#293](https://github.com/gridap/Gridap.jl/pull/293).
- New `FESpace` constructor that takes an instance of `ReferenceFE`. Since PR [#294](https://github.com/gridap/Gridap.jl/pull/294).
- New `FESpace` constructor that takes an instance of `Conformity`. Since PR [#311](https://github.com/gridap/Gridap.jl/pull/311).
- New `CDLagrangianRefFE` struct, that provides a Lagrangian reference FE with different conformity per direction. Since PR [#299](https://github.com/gridap/Gridap.jl/pull/299).
- New `FESpace` method that takes a model and a `RefFE`. Since PR [#299](https://github.com/gridap/Gridap.jl/pull/299).
- Possibility to have 0 order in `DISC` directions of a `CDLagrangianRefFE`. Since PR [#308](https://github.com/gridap/Gridap.jl/pull/308).
- Added setindex! method for Reindexed. Since PR [#309](https://github.com/gridap/Gridap.jl/pull/309).

### Changed

- Changed the interfaces of `ReferenceFE` and `NodalReferenceFE` in relation of DOF ownership. Now function `get_face_own_dofs` and related ones are parametrized by a `Conformity` object. Since PR [#311](https://github.com/gridap/Gridap.jl/pull/311).
- The constructors `GenericRefFE`, `GenericNodalCartesianRefFE`, and `compute_conforming_cell_dofs` take an extra argument of type `Conformity`. Since PR [#311](https://github.com/gridap/Gridap.jl/pull/311).
- Renamed `PDiscRefFE` -> `DiscRefFE` struct keeping the name for constructor. Since PR [#293](https://github.com/gridap/Gridap.jl/pull/293).
- One of the `GradConformingFESpace` methods now more general `ConformingFESpace`. Since PR [#293](https://github.com/gridap/Gridap.jl/pull/293).
- `DivConformingFESpace` and `CurlConformingFESpace` constructors eliminated. Since PR [#293](https://github.com/gridap/Gridap.jl/pull/293).
- Extend table to support arbitrary vector types. Since PR [#310](https://github.com/gridap/Gridap.jl/pull/310).

### Fixed

- Construction of `VectorValue`, `TensorValue`, et al. objects from non-homogeneous arguments. This solves some problems associated with automatic differentiation. Since PR [#298](https://github.com/gridap/Gridap.jl/pull/298).
- `CDLagrangianRefFE` node ordering. Since PR [#305](https://github.com/gridap/Gridap.jl/pull/305).

## [0.11.2] - 2020-06-22

### Added

- Method `solve!(x,ls,op::AffineOperator,cache::Nothing,newmatrix)`. Since PR [#288](https://github.com/gridap/Gridap.jl/pull/288).

### Fixed

- Bug related with `WriteVTK` version 1.7. Fixed via PR [#287](https://github.com/gridap/Gridap.jl/pull/287).
- Bug in outer constructor of Table{...} for input arrays of abstract type. Fixed via PR [#285](https://github.com/gridap/Gridap.jl/pull/285).

## [0.11.1] - 2020-06-19

### Fixed

- Bug in the handling of caches in `NLSolver`. Fixed via PR [#283](https://github.com/gridap/Gridap.jl/pull/283).
- Bug that showed up when interpolating a FE function defined on an
  `ExtendedFESpace` onto a non-extended `FESpace`. Fixed via PR [#282](https://github.com/gridap/Gridap.jl/pull/282).

## [0.11.0] - 2020-06-16

### Added

- Operator `⊙` (\odot) as an alias of `inner`. Since PR [#239](https://github.com/gridap/Gridap.jl/pull/239).
- Operator `⊗` (\otimes) as an alias of `outer`. Since PR [#239](https://github.com/gridap/Gridap.jl/pull/239).
- Support for (symmetric) 4th order tensors. Since PR [#239](https://github.com/gridap/Gridap.jl/pull/239).
- Optimizations for symmetric 2nd order tensors. Since PR [#239](https://github.com/gridap/Gridap.jl/pull/239).
- Methods for `cross` function (aka `×` (\times)) to operate with `VectorValues`. Since PR [#280](https://github.com/gridap/Gridap.jl/pull/280).
- Interpolation is now supported also for multifield spaces. Since PR [#279](https://github.com/gridap/Gridap.jl/pull/279).

### Changed

- Major refactoring in the module `Gridap.TensorValues`.
  Since PR [#239](https://github.com/gridap/Gridap.jl/pull/239).
  **The following changes are likely to affect all users:**
  - The operator `*` is not allowed for expressing the dot product anymore. Use `LinearAlgebra.dot`
    function aka `⋅` (\cdot).
  - The syntax `∇*u` is not allowed anymore. Use `∇⋅u` instead.
  - Gridap re-exports `dot`, `⋅`, and other names from LinearAlgebra that are used
    often in Gridap code.
  - Function `n_components` is renamed to `num_components`.
- The `SingleFieldFESpace` interface has changed. The function `gather_free_and_dirichlet_values!`
  has been added as mandatory for all FE space implementations and the old function `gather_free_and_dirichlet_values`
  is now optional. Since PR [#279](https://github.com/gridap/Gridap.jl/pull/279).

## [0.10.4] - 2020-06-8

### Added

- Functions `create_vtk_file` and `createvtk`. Since PR [#273](https://github.com/gridap/Gridap.jl/pull/273).

## [0.10.3] - 2020-05-29

### Added

- Function `print_op_tree` to visualize lazy operation trees. Since PR [#270](https://github.com/gridap/Gridap.jl/pull/270).
- Exported `apply` and `reindex` from `Gridap` top level. Since PR [#270](https://github.com/gridap/Gridap.jl/pull/270).
- Extended support of `CartesianDiscreteModel` to models with periodic boundary conditions.
  PR [#266](https://github.com/gridap/Gridap.jl/pull/266).

### Deprecated

- Optional argument `map` for CartesianDescriptor converted to a key-word argument. Since PR [#266](https://github.com/gridap/Gridap.jl/pull/266).

### Fixed

- Fixed some methods of the `sparsecsr` generic function. Since PR [#262](https://github.com/gridap/Gridap.jl/pull/262).
- Fixed BUG in `findnz` function for `SparseMatrixCSR`. Since PR [#264](https://github.com/gridap/Gridap.jl/pull/264).
- Fixed `restrict(::AbstractArray,::TriangulationPortion)` for portions of triangulations extending `BoundaryTriangulation`. Since PR [#267](https://github.com/gridap/Gridap.jl/pull/267).

## [0.10.2] - 2020-05-21

### Added

- New key-word arguments `zeromean_trian` and `zeromean_quad` in the `FESpace` constructor. Since
  PR [#257](https://github.com/gridap/Gridap.jl/pull/257).
- New method `reindex(::Triangulation,indices)`. Since
  PR [#257](https://github.com/gridap/Gridap.jl/pull/257).
- New functions `get_face_to_face(::BoundaryTriangulation)` and `get_cell_around(::BoundaryTriangulation)`. Since
  PR [#256](https://github.com/gridap/Gridap.jl/pull/256).

## [0.10.1] - 2020-05-19

### Fixed

- Added missing implementation of `simplexify(SEGMENT)` and `simplexify(VERTEX)`. Since PR [#252](https://github.com/gridap/Gridap.jl/pull/252).

## [0.10.0] - 2020-05-14

### Added

- Extended support of `TriangulationPortion` to boundary and skeleton triangulations. Since PR [#249](https://github.com/gridap/Gridap.jl/pull/249).
- Added `FESpaceWithLinearConstraints`. Since PR [#247](https://github.com/gridap/Gridap.jl/pull/247).
- Added inner constructor to `CartesianDiscreteModel` allowing to build a model that represents a subgrid of
  a larger grid. Since PR [#245](https://github.com/gridap/Gridap.jl/pull/245).

### Changed

- The part associated with the imposition of constraints in the `FESpace` interface has changed slightly. Since PR [#247](https://github.com/gridap/Gridap.jl/pull/247).
- Simplified the signature of `zero_free_values(::FESpace)`. Since PR [#249](https://github.com/gridap/Gridap.jl/pull/249).
- Simplified the signature of `zero_initial_guess(op::NonlinearOperator)`. Since PR [#249](https://github.com/gridap/Gridap.jl/pull/249).
- Major refactoring in the `Assembler` interface.
  **Important change:** Now, assembly-related functions take the data returned by functions like
  `collect_cell_matrix` as it is. Example: the old user code `assemble_matrix(assembler,collect_cell_matrix(du,dv,terms)...)`
  now is written simply as `assemble_matrix(assembler,collect_cell_matrix(du,dv,terms))`, i.e., the unpack of the last argument is not
  used anymore. In addition, with the new assembler interface, it is possible to customize the assembly process
  via a so-called `AssemblerStrategy` object. Since PR [#249](https://github.com/gridap/Gridap.jl/pull/249).
- Change the types of the sizes and partition fields of CartesianDescriptor to tuples instead of points.
  Since PR [#246](https://github.com/gridap/Gridap.jl/pull/246).

## [0.9.2] - 2020-04-26

### Added

- Automatic differentiation of manufactured solutions. Since PR [#236](https://github.com/gridap/Gridap.jl/pull/236).

## [0.9.1] - 2020-04-20

### Added

- Function `cell_measure`. Since PR [#234](https://github.com/gridap/Gridap.jl/pull/234).

### Fixed

- Several bugs associated with `ExtendedFESpace`. In particular, we have fixed a bug that showed up when combining `ZeroMeanFESpace` and `ExtendedFESpace`. Since PR [#234](https://github.com/gridap/Gridap.jl/pull/234).

## [0.9.0] - 2020-04-18

### Added

- Function `HomogeneousTrialFESpace`. Since PR [#226](https://github.com/gridap/Gridap.jl/pull/226).
- Function `lazy_append` in order to lazily append two objects (implemented for `AbstractVector`, `Triangulation`, and `CellQuadrature`). Since PR [#220](https://github.com/gridap/Gridap.jl/pull/220).
- Support for FE spaces with DOFs defined in the physical space. Since PR [#216](https://github.com/gridap/Gridap.jl/pull/216) and [#218](https://github.com/gridap/Gridap.jl/pull/218).

### Changed

- Replaced `non_linear` -> `nonlinear` and `NonLinear` -> `Nonlinear`. Since PR [#223](https://github.com/gridap/Gridap.jl/pull/223).
- The `FESpace` interface has slightly changed, mainly the return type of functions `get_cell_basis` and `get_cell_dof_basis.`. Since PR [#216](https://github.com/gridap/Gridap.jl/pull/216) and [#218](https://github.com/gridap/Gridap.jl/pull/218).

### Fixed

- Bug that showed up in multi-field computations when some field had no contribution to the rhs vector. Since [#229](https://github.com/gridap/Gridap.jl/pull/229).
- Bug in gradient operator in the void part of `ExtendedFESpace` objects. Since PR [#219](https://github.com/gridap/Gridap.jl/pull/219).
- Bug in jumps of quantities restricted to `InterfaceTriangulation` objects. Since PR [#215](https://github.com/gridap/Gridap.jl/pull/215).

## [0.8.0] - 2020-03-17

### Added

- Support for surface-coupled multi-physics. See [`SurfaceCouplingTests.jl`](https://github.com/gridap/Gridap.jl/blob/master/test/GridapTests/SurfaceCouplingTests.jl) for further details. Since PR [#209](https://github.com/gridap/Gridap.jl/pull/209).
- Support for constitutive laws with state / historical variables. See [`IsotropicDamageTests.jl`](https://github.com/gridap/Gridap.jl/blob/master/test/GridapTests/IsotropicDamageTests.jl) for further details. Since PR [#208](https://github.com/gridap/Gridap.jl/pull/208).
- Curl-conforming reference FE `NedelecRefFE` and corresponding FE space constructor since PR [#199](https://github.com/gridap/Gridap.jl/pull/199).
- New constructors `AffineFETermFromCellMatVec` and `FETermFromCellJacRes` that provides full control in the definition of cell matrices and vectors. Since PR [#191](https://github.com/gridap/Gridap.jl/pull/191).
- Support for simultaneous integration of matrices and vectors. Since PR [#191](https://github.com/gridap/Gridap.jl/pull/191).

### Changed

- Renaming NonLinear to Nonlinear since it is one word and it is not consistent with style
- The definition of interpolation order in Raviart-Thomas and Nédélec reference FEs has changed. Now, the divergence of functions in the Raviart-Thomas space of order `k` belongs to `P_k` or `Q_k` depending on the underlying polytope. Idem for Nédelec, but using the curl instead of the divergence. Since PR [#212](https://github.com/gridap/Gridap.jl/pull/212).

- The order in which test and trial spaces are written in the code has changed and also the other in the arguments of functions defining bi-linear and linear forms, and weak residuals and Jacobians. **This affects everybody that is using Gridap, even the most basic users**. Now, we write the trial space before the test one in all methods taking two spaces in their arguments. E.g., we have changed `AffineFEOperator(V,U,terms...)` to `AffineFEOperator(U,V,terms...)`, where `U` is the trial and `V` is the test space. For functions defining weak forms, now we have: The new signatures for bi-linear and a linear forms are `a(u,v)`, `l(v)`, where `u` is a trial function and `v` is a test one. For weak Jacobians and residuals `jac(u,du,v)` and `res(u,v)`, where `u` is the (trial) function in which we evaluate these quantities, `du` is the direction in which we evaluate the Jacobian and `v` is a test function. Since PR [#195](https://github.com/gridap/Gridap.jl/pull/195) and PR [#197](https://github.com/gridap/Gridap.jl/pull/197).

- The part related with the application of constraints in the `FESpace` interface has changed. Since PR [#191](https://github.com/gridap/Gridap.jl/pull/191).

### Fixed

- Bug in 1d Cartesian grids. Since PR [#192](https://github.com/gridap/Gridap.jl/pull/192).

## [0.7.1] - 2020-02-18

### Added

- New `DirichletFESpace` that can be used to compute matrices and vectors associated with the Dirichlet DOFs. Since commit [972afcc](https://github.com/gridap/Gridap.jl/commit/972afcc6dd8e024a7daeebd160a9dabe44ff5921)

## [0.7.0] - 2020-02-13

This version is a major refactoring of the project which is not summarized here for the sake of brevity. Most of the functionality of v0.6.0 is available in v0.7.0, but with a possibly slightly different API. See [here](https://github.com/gridap/Tutorials/compare/v0.6.0...v0.7.0) the changes in the sources of the Gridap Tutorials between versions 0.6.0 and 0.7.0 to effectively see the major changes in the API.

## [0.6.0] - 2020-01-24

### Added

- New `GenericRefFE`. Since commit [876ef1e](https://github.com/gridap/Gridap.jl/commit/c3c9010177432b8f07aaecf4a0baa4b93876ef1e)
- New `NedelecRefFE` constructor that generates Nedelec FEs of arbitrary order in 2D and 3D on hex. Since commit [876ef1e](https://github.com/gridap/Gridap.jl/commit/c3c9010177432b8f07aaecf4a0baa4b93876ef1e)
- New keyword argument `map` in the constructor of `CartesianModel`, which allows one to transform the original domain, by default [0,1]^d to a new domain through a homeomorphic map. Since commit [30cc4bc](https://github.com/gridap/Gridap.jl/commit/1c51b18f7e21c4915c0b379585dc5d98130cc4bc)
- New keyword argument `map` in the constructor of `CartesianGrid` and a new `map` attribute in this structure, since commit [30cc4bc](https://github.com/gridap/Gridap.jl/commit/1c51b18f7e21c4915c0b379585dc5d98130cc4bc)
- `CartesianGridPoints` has new attribute `map` since commit [30cc4bc](https://github.com/gridap/Gridap.jl/commit/1c51b18f7e21c4915c0b379585dc5d98130cc4bc)
- Added [`SparseMatricesCSR`](https://github.com/gridap/SparseMatricesCSR.jl) support to `SparseMatrixAssembler` and `MultiSparseMatrixAssembler` in [PR #118](https://github.com/gridap/Gridap.jl/pull/118#).

### Changed

- The `RaviartThomasRefFE` has now been replaced by `GenericRefFE`, and the constructor for Raviart-Thomas FEs is called `RTRefFE`. Since commit [876ef1e](https://github.com/gridap/Gridap.jl/commit/c3c9010177432b8f07aaecf4a0baa4b93876ef1e)
- The default map in the `CartesianModel` constructor is [0,1]^d instead of [-1,1]^d, since commit [30cc4bc](https://github.com/gridap/Gridap.jl/commit/1c51b18f7e21c4915c0b379585dc5d98130cc4bc)
- `CartesianGrid` has attribute `map` since commit [30cc4bc](https://github.com/gridap/Gridap.jl/commit/1c51b18f7e21c4915c0b379585dc5d98130cc4bc)

## [0.5.2] - 2019-10-22

### Fixed

- Incompatibility problem with `TensorValues` version 0.3.5. Via commit [3c0682a](https://github.com/gridap/Gridap.jl/commit/3c0682a84250e17086457ca3a90a49d9bce133d0).

## [0.5.1] - 2019-10-03

### Added

- Pretty printing for the types most exposed to users. Since PR [#109](https://github.com/gridap/Gridap.jl/pull/109).

### Fixed

- Bug related to `ZeroMeanFESpace`. Via PR [#111](https://github.com/gridap/Gridap.jl/pull/111).

## [0.5.0] - 2019-09-27

### Added

- Added a high level constructor, namely `FESpace`, to create different types of FE spaces. See issue [#100](https://github.com/gridap/Gridap.jl/issues/100) for more details. Since PR [#102](https://github.com/gridap/Gridap.jl/pull/102).
- Added `ZeroMeanFESpace` to construct FE spaces whose functions have zero mean value. Since PR [#102](https://github.com/gridap/Gridap.jl/pull/102).
- Added Hdiv FE space using Raviart-Thomas reference FEs in [34bfa34](https://github.com/gridap/Gridap.jl/commit/34bfa344efd1bc6a5d3c5993d9639259ed21671a)
- Added the corresponding DOF basis for Raviart-Thomas reference FEs for interpolation of fields [60b9021](https://github.com/gridap/Gridap.jl/commit/60b9021b6d4b5e66a9ec4fe2067aa8278f8ccb52)
- Added an arbitrary order div-conforming Raviart-Thomas reference FE of arbitrary order on quads in commit
  [60b9021](https://github.com/gridap/Gridap.jl/commit/60b9021b6d4b5e66a9ec4fe2067aa8278f8ccb52)
- Now, the `tags` argument is optional when constructing `SkeletonTriangulation` and `BoundaryTriangulation` objects from a `DiscreteModel`. Since commit [e6424a3](https://github.com/gridap/Gridap.jl/commit/e6424a304feb38547241e86de07a821e26344a7e).
- Added `mean` operator for quantities restricted to a `SkeletonTriangulation`. Since commit [83798b4](https://github.com/gridap/Gridap.jl/commit/83798b4f38aaf482b968ffd0359eb75c79a21385).
- Extended `NormalVector` to `SkeletonTriangulations`. Since commit [5fb8487](https://github.com/gridap/Gridap.jl/commit/5fb84871128c4388559cc5052d9ff00f0be19462).
- Now, `TrialFESpaces` can be constructed from values instead of functions if the corresponding Dirichlet conditions are constant. Since commit [bae237e](https://github.com/gridap/Gridap.jl/commit/bae237e881db6569622f3559f82bcc3999560526).
- Added the possibility of adding new tags to a `FaceLabels` object via the function `add_tag_from_tags!` and using it to construct FE spaces. Since commit [e9dfac4](https://github.com/gridap/Gridap.jl/commit/e9dfac4489047c0b7e1c62507f4335e9fc76dfd8).
- Added `BackslashSolver` to facilitate the usage in Gridap of the build-in Julia backslash linear solver. Since commit [8e3a9b7](https://github.com/gridap/Gridap.jl/commit/8e3a9b71c64b032c5a572a7ef696f4cbf875190b).
- Added `NLSolver` to facilitate the usage in Gridap of the non-linear solvers available in the official Julia package `NLsolve`. Introduced in commit [e5a933f](https://github.com/gridap/Gridap.jl/commit/e5a933f3093faea221a50bdd796d7f02113ed52c) as `JuliaNLSolver`. Renamed to `NLSolver` in PR [#108](https://github.com/gridap/Gridap.jl/pull/108).

### Changed

- The Signature of `solve!` for `NumericalSetup` objects. The argument for the system matrix has been removed. The information about the matrix is already in the `NumericalSetup` object. Since commit [ac212d3](https://github.com/gridap/Gridap.jl/commit/ac212d30205700a919a37f9abf9dac6cbde03e38).
- The signature of `solve!(::FEFunction,::FESolver,::FEOperator)`. Before it was used as `cache = solve!(uh,solver,op)`, now it is used as `uh, cache = solve!(uh,solver,op)`. Since PR [#102](https://github.com/gridap/Gridap.jl/pull/102).
- Previous ConformingFESpace constructor is H1ConformingFESpace since [34bfa34](https://github.com/gridap/Gridap.jl/commit/34bfa344efd1bc6a5d3c5993d9639259ed21671a)

### Deprecated

- `JuliaNLSolver`. Renamed to `NLSolver`. Since PR [#108](https://github.com/gridap/Gridap.jl/pull/108).
- Key-word argument `order` in `CellQuadrature` constructor. Renamed to `degree`. Since PR [#108](https://github.com/gridap/Gridap.jl/pull/108).

### Fixed

- Bug in `@law` macro for more than one `FEBasis` arguments. Solved via PR [#104](https://github.com/gridap/Gridap.jl/pull/104).
- Bug in `NonlinearFEOperator` constructor with default assembler in multi-field computations. Solved via PR [#104](https://github.com/gridap/Gridap.jl/pull/104).
- Bug in `NormalVector` for non-Cartesian grids. Solved via PR [#98](https://github.com/gridap/Gridap.jl/pull/98).

## [0.4.0] - 2019-09-07

### Added

- Added support to high order simplicial Lagrangian finite elements. Since commit [cbefe9b](https://github.com/gridap/Gridap.jl/commit/cbefe9bbea83d00e7f6ccbef50396ddc7dc49b80).
- Now the built-in simplicial grids are oriented. Since commit [cbefe9b](https://github.com/gridap/Gridap.jl/commit/cbefe9bbea83d00e7f6ccbef50396ddc7dc49b80).
- Added binary operations between `FEFunction` and `Number`, and `FEBasis` and `Number`. Since PR [#88](https://github.com/gridap/Gridap.jl/pull/88).
- Added `PDiscRefFE`, `DiscFESpace`, and `ConstrainedFESpace`. Since PR [#88](https://github.com/gridap/Gridap.jl/pull/88).
- Now its possible to pass a `CellNumber` or an `Array` of numbers into a constitutive law. Useful to identify which is the material of the current Gauss point in multi-material problems. Since commit [62cb2c3](https://github.com/gridap/Gridap.jl/commit/62cb2c354e2a09c556324a4fe9861329989299f4).
- `LinearFESolver` is now optional for solving a `LinearFEOperator`. Since commit [5c1caa8](https://github.com/gridap/Gridap.jl/commit/5c1caa8c92b260db72f5902e778ec5c0eb88728b).
- `Assembler` is now optional to build `FEOperator` objects. Since commit [b1bf517](https://github.com/gridap/Gridap.jl/commit/b1bf5172955b940f6b3c9d027bd4a839c6486199).
- Binary operations between `Function` and `FEFunction`. Since commit [a7f22f5](https://github.com/gridap/Gridap.jl/commit/a7f22f5ac1f45d9e8f53906472257aa582726e87).
- Extended constructions of `CLagrangianFESpace` and `DLagrangianFESpace`. `diritags` and `dirimasks` are now optional. `diritags` can now be also a vector of `String`. Since commit [776b402](https://github.com/gridap/Gridap.jl/commit/776b40238365f145037fc5e490600bf5b45434ef).
- Added `div`, `curl`, and `trace` operators. Since commit [5a0f322](https://github.com/gridap/Gridap.jl/commit/5a0f322c5b938f12e26e9c0a7c9361aa649e014f).
- Macro `@law` to facilitate the definition of constitutive laws. Since commit [30b67f2](https://github.com/gridap/Gridap.jl/commit/30b67f29009b872944be94486dc4a1b0134a0a60).
- Definition of linear forms `b(v) = inner(v, f)` directly from a function `f`. Since commit [bb42847](https://github.com/gridap/Gridap.jl/commit/bb42847c702a99b9b5f2c2d922fbe4c95b23f646)
- Serialization and de-serialization of `DiscreteModel` objects into and from `json` format. Since PR [#76](https://github.com/gridap/Gridap.jl/pull/76).
- Support for boundary integration (e.g., Neumann BCs) for multi-field computations. Since PR [#75](https://github.com/gridap/Gridap.jl/pull/75).

### Changed

- Signature of `LagrangianRefFE` constructor. Since commit [529c764](https://github.com/gridap/Gridap.jl/commit/529c7646a531db6910a00f04a925dadec3a50b7c).

### Fixed

- Bug in `LinearFETerm` for multi-field computations. Fixed via commit [2b957d1](https://github.com/gridap/Gridap.jl/commit/2b957d1b3a9a9a4396075801d8c837f6aff921c8).
- Bug in `MultiCellArray` constructor. Fixed via commit [bbc3b1c](https://github.com/gridap/Gridap.jl/commit/bbc3b1c91752f8efa978731cb90c6198dc0e5227).
- Bug in binary operations between FEFunction and FEBasis. Fixed via commit [aa49689](https://github.com/gridap/Gridap.jl/commit/aa49689be2a8dc14e052a6409c8348f492b52b3e).

## [0.3.0] - 2019-08-06

### Added

- `CurlGradMonomialBasis` spanning the polynomial space needed for RT elements on n-cubes.
- `CLagrangianFESpace` and `DLagrangianFESpace` types providing an efficient implementation for continuous and discontinuous Lagrangian FE spaces respectively. In contrast to `ConfirmingFESpace`, the new types allow to select which are the components that are actually prescribed on the Dirichlet boundary. Since PR [#64](https://github.com/gridap/Gridap.jl/pull/64).
- `simplexify` function to convert `Grid` and `DiscreteModel` objects made of n-cubes to the corresponding counterparts made of n-simplices. Since PR [#62](https://github.com/gridap/Gridap.jl/pull/62).
- Duffy transformation based integration for n-simplices of arbitrary dimension. Since PR [#61](https://github.com/gridap/Gridap.jl/pull/61).
- `NormalVector` to construct the outward normal vector to a given `BoundaryTriangulation`. Since PR [#60](https://github.com/gridap/Gridap.jl/pull/60).
- Support for tensor-valued FE computations. Since PR [#57](https://github.com/gridap/Gridap.jl/pull/57).
- Support for integration on the skeleton of the mesh. This includes `SkeletonTriangulation`, an integration mesh for the skeleton, `restrict` function is extended to restrict to the skeleton, `jump` function to compute jumps of `CellFields` and `CellBasis` restricted to the skeleton, extension of `FETerms` to allow integration on the skeleton. See PR [#47](https://github.com/gridap/Gridap.jl/pull/47)
- Support for Robin boundary conditions. Since commit [946054a](https://github.com/gridap/Gridap.jl/commit/946054a028e658afa87c7a7c71e973957a2c4877)
- Support for Neumann boundary conditions. Since commit [4dcd16f](https://github.com/gridap/Gridap.jl/commit/4dcd16fbff9edf66fb66efb748ef01901c20a4aa)
- `FETerm` and `AffineFETerm` abstract types and several concrete implementations. They allow to deal with problems whose weak form has terms integrated over different geometrical entities. `NonlinearFEOperator` and `LinearFEOperator` can be constructed using several terms. Since commit [0f99234](https://github.com/gridap/Gridap.jl/commit/0f99234156dd0174485ca83431de76aa3825584a)
- Extended `Assembler` and `MultiAssembler` to deal with several terms. See issue [#42](https://github.com/gridap/Gridap.jl/issues/42) and PR [#43](https://github.com/gridap/Gridap.jl/pull/43).
- `IdentityCellNumber`, an indexable cell number that simply returns the given index. Also efficient implementation of `reindex` for this type (i.e. do nothing). Available since commit [b6b4c32](https://github.com/gridap/Gridap.jl/commit/b6b4c32c8c4b826a41ba64c770ac8a1c394e16f0)
- Function `restrict` for restricting `CellField` and `CellBasis` objects to surfaces. Available since commit [e981f3c](https://github.com/gridap/Gridap.jl/commit/e981f3c221f3624cfc6764efa47f22652fc22b4f)
- `BoundaryTriangulation` an integration mesh used to integrate `CellField` and `CellBasis` objects restricted on a surface. Available since commit [e981f3c](https://github.com/gridap/Gridap.jl/commit/e981f3c221f3624cfc6764efa47f22652fc22b4f)
- `NonIterableCellMap`, a cell map that has iteration intentionally disabled. Available since commit [956a537](https://github.com/gridap/Gridap.jl/commit/956a5374db6c3b9546e85e0d4d49ae0560057565).
- `CompressedCellValue`, `CompressedCellArray`, and `CompressedCellMap`, as well as efficient versions of `apply`, `evaluate`, and `reindex` for these types. See PR [#41](https://github.com/gridap/Gridap.jl/pull/41) for more details.
- `NEWS.md` file (a changelog file)

### Changed

- Domains are now in [0,1] instead of [-1,1]. Quadratures and nodes arrays modified accordingly. Since commit [268dfe1](https://github.com/gridap/Gridap.jl/commit/268dfe12ef7d736fcd9ad0b9b256740aaf15b2e7).
- Changed the signature of `assemble`, `apply_constraints`, `apply_constraints_rows`, and `apply_constraints_cols` to support FE assembly of several terms, which are integrated in different domains. The old API of `asseble` is still functional, but not for the `apply_constraints` et al. Since PR [#43](https://github.com/gridap/Gridap.jl/pull/43). Further changed in commit [a335aed](https://github.com/gridap/Gridap.jl/commit/a335aede65c92a1f61f0ff0dbb0fb44cc20cf906).

### Fixed

- Bug in generation of the cellwise local to global DOF map for high order interpolations. Fixed via PR [#56](https://github.com/gridap/Gridap.jl/pull/56).
- Bug in numerical integration. There was a bug for computations where the number of cell DOFs was different from the number of integration points. Fixed via commit [0b3d4bf](https://github.com/gridap/Gridap.jl/commit/0b3d4bfadea48707c748fca0de65a51a598b6ca6)

## [0.2.0] - 2019-06-29

A changelog is not maintained for this version.

This version introduces the core finite element machinery for linear and non-linear problems,
single field and multi-field problems with terms integrated over the interior of the computational domain.

## [0.1.0] - 2019-05-20

A changelog is not maintained for this version.

This version is non functional. It is just a tag for registering the package.
