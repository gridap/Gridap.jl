# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added Hdiv FE space using Raviart-Thomas reference FEs in [34bfa34](https://github.com/gridap/Gridap.jl/commit/34bfa344efd1bc6a5d3c5993d9639259ed21671a)
- Added the corresponding DOF basis for Raviart-Thomas reference FEs for interpolation of fields [60b9021](https://github.com/gridap/Gridap.jl/commit/60b9021b6d4b5e66a9ec4fe2067aa8278f8ccb52)
- Added an arbitrary order div-conforming Raviart-Thomas reference FE of arbitrary order on quads in commit
[60b9021](https://github.com/gridap/Gridap.jl/commit/60b9021b6d4b5e66a9ec4fe2067aa8278f8ccb52)
- Now, the `tags` argument is optional when constucting `SkeletonTriangulation` and `BoundaryTriangulation` objects from a `DiscreteModel`. Since commit [e6424a3](https://github.com/gridap/Gridap.jl/commit/e6424a304feb38547241e86de07a821e26344a7e).
- Added `mean` operator for quantities restricted to a `SkeletonTriangulation`. Since commit [83798b4](https://github.com/gridap/Gridap.jl/commit/83798b4f38aaf482b968ffd0359eb75c79a21385).  
- Extended `NormalVector` to `SkeletonTriangulations`. Since commit [5fb8487](https://github.com/gridap/Gridap.jl/commit/5fb84871128c4388559cc5052d9ff00f0be19462). 
- Now, `TrialFESpaces` can be constructed from values instead of functions if the corresponding Dirichlet conditions are constant. Since commit [bae237e](https://github.com/gridap/Gridap.jl/commit/bae237e881db6569622f3559f82bcc3999560526).
- Added the possibility of adding new tags to a `FaceLabels` object via the function `add_tag_from_tags!` and using it to construct FE spaces. Since commit [e9dfac4](https://github.com/gridap/Gridap.jl/commit/e9dfac4489047c0b7e1c62507f4335e9fc76dfd8).
- Added `BackslashSolver` to facilitate the usage in Gridap of the build-in Julia backslash linear solver. Since commit [8e3a9b7](https://github.com/gridap/Gridap.jl/commit/8e3a9b71c64b032c5a572a7ef696f4cbf875190b).
- Added `JuliaNLSolver` to facilitat the usage in Gridap of the non-linear solvers available in the official Julia package `NLsolve`. Since commit [e5a933f](https://github.com/gridap/Gridap.jl/commit/e5a933f3093faea221a50bdd796d7f02113ed52c).

### Changed

- Previous ConformingFESpace constructor is H1ConformingFESpace since [34bfa34](https://github.com/gridap/Gridap.jl/commit/34bfa344efd1bc6a5d3c5993d9639259ed21671a)

### Removed
### Deprecated
### Fixed

## [0.4.0] - 2019-09-07

### Added

- Added support to high order simplicial Lagrangian finite elements. Since commit [cbefe9b](https://github.com/gridap/Gridap.jl/commit/cbefe9bbea83d00e7f6ccbef50396ddc7dc49b80).
- Now the built-in simplicial grids are oriented. Since commit [cbefe9b](https://github.com/gridap/Gridap.jl/commit/cbefe9bbea83d00e7f6ccbef50396ddc7dc49b80).
- Added binary operations between `FEFuntion` and `Number`, and `FEBasis` and `Number`.  Since PR [#88](https://github.com/gridap/Gridap.jl/pull/88).
- Added `PDiscRefFE`, `DiscFESpace`, and `ConstrainedFESpace`. Since PR [#88](https://github.com/gridap/Gridap.jl/pull/88).
- Now its possible to pass a `CellNumer` or an `Array` of numbers into a constitutive law. Usefull to identify which is the material of the current Gauss point in multi-material problems. Since commit [62cb2c3](https://github.com/gridap/Gridap.jl/commit/62cb2c354e2a09c556324a4fe9861329989299f4).
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
