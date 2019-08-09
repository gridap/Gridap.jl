# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Definition of linear forms `b(v) = inner(v, f)` directly from a function `f`. Since commit [bb42847](https://github.com/gridap/Gridap.jl/commit/bb42847c702a99b9b5f2c2d922fbe4c95b23f646)
- Serialization and de-serialization of `DiscreteModel` objects into and from `json` format. Since PR [#76](https://github.com/gridap/Gridap.jl/pull/76).
- Support for boundary integration (e.g., Neumann BCs) for multi-field computations. Since PR [#75](https://github.com/gridap/Gridap.jl/pull/75).

### Changed
### Removed
### Deprecated
### Fixed

## [0.3.0] - 2019-08-06
### Added
- `CurlGradMonomialBasis` spanning the polynomial space needed for RT elements on n-cubes
- `CLagrangianFESpace` and `DLagrangianFESpace` types providing an efficient implementation for continuous and discontinuous Lagrangian FE spaces respectivelly. In contrast to `ConfirmingFESpace`, the new types allow to select which are the components that are actually prescribed on the Dirichlet boundary. Since PR [#64](https://github.com/gridap/Gridap.jl/pull/64).
- `simplexify` funciton to convert `Grid` and `DiscreteModel` objects made of n-cubes to the corresponding counterparts made of n-simplices. Since PR [#62](https://github.com/gridap/Gridap.jl/pull/62).
- Duffy transformation based integration for n-simplices of arbitrary dimension. Since PR [#61](https://github.com/gridap/Gridap.jl/pull/61).
- `NormalVector` to construct the outward normal vector to a given `BoundaryTriangulation`. Since PR [#60](https://github.com/gridap/Gridap.jl/pull/60).
- Support for tensor-valued FE computations. Since PR [#57](https://github.com/gridap/Gridap.jl/pull/57).
- Support for integration on the skeleton of the mesh. This includes `SkeletonTriangulation`, an integration mesh for the skeleton, `restrict` function is extended to restrict to the skeleton, `jump` function to compute jumps of `CellFields` and `CellBasis` restricted to the skeleton, extension of `FETerms` to allow integration on the skeleton. See PR [#47](https://github.com/gridap/Gridap.jl/pull/47)
- Support for Robin boundary conditions. Since commit [946054a](https://github.com/gridap/Gridap.jl/commit/946054a028e658afa87c7a7c71e973957a2c4877)
- Support for Neumann boundary conditions. Since commit [4dcd16f](https://github.com/gridap/Gridap.jl/commit/4dcd16fbff9edf66fb66efb748ef01901c20a4aa)
- `FETerm` and `AffineFETerm` abstract types and several concrete implementations. They allow to deal with problems whose weak form has terms integrated over different geometrical entities. `NonLinearFEOperator` and `LinearFEOperator` can be constructed using several terms. Since commit [0f99234](https://github.com/gridap/Gridap.jl/commit/0f99234156dd0174485ca83431de76aa3825584a)
- Extended `Assembler` and `MultiAssembler` to deal with several terms. See issue [#42](https://github.com/gridap/Gridap.jl/issues/42) and PR [#43](https://github.com/gridap/Gridap.jl/pull/43).
- `IdentityCellNumber`, an indexable cell number that simply returns the given index. Also efficient implementation of `reindex` for this type (i.e. do nothing). Available since commit [b6b4c32](https://github.com/gridap/Gridap.jl/commit/b6b4c32c8c4b826a41ba64c770ac8a1c394e16f0)
- Function `restrict` for restricting `CellField` and `CellBasis` objects to surfaces. Available since commit [e981f3c](https://github.com/gridap/Gridap.jl/commit/e981f3c221f3624cfc6764efa47f22652fc22b4f)
- `BoundaryTriangulation` an integration mesh used to integrate `CellField` and `CellBasis` objects restricted on a surface. Available since commit [e981f3c](https://github.com/gridap/Gridap.jl/commit/e981f3c221f3624cfc6764efa47f22652fc22b4f)
- `NonIterableCellMap`, a cell map that has iteration intentionally disabled. Available since commit [956a537](https://github.com/gridap/Gridap.jl/commit/956a5374db6c3b9546e85e0d4d49ae0560057565).
- `CompressedCellValue`, `CompressedCellArray`, and `CompressedCellMap`, as well as efficient versions of `apply`, `evaluate`, and `reindex` for these types. See PR [#41](https://github.com/gridap/Gridap.jl/pull/41) for more details.
- `NEWS.md` file (a changelog file)

### Changed
- Domains are now in [0,1] instead of [-1,1]. Quadratures and nodes arrays modified accordingly. Since commit [268dfe1](https://github.com/gridap/Gridap.jl/commit/268dfe12ef7d736fcd9ad0b9b256740aaf15b2e7).
- Changed the signature of `assemble`, `apply_constraints`, `apply_constraints_rows`, and `apply_constraints_cols` to support FE assembly of several terms, which  are integrated in different domains. The old API of `asseble` is still functional, but not for the `apply_constraints` et al. Since PR [#43](https://github.com/gridap/Gridap.jl/pull/43). Further changed in commit   [a335aed](https://github.com/gridap/Gridap.jl/commit/a335aede65c92a1f61f0ff0dbb0fb44cc20cf906).

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

