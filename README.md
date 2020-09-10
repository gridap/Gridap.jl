# <img src="https://github.com/gridap/Gridap.jl/blob/master/images/color-text.png" width="250" title="Gridap logo">


| **Documentation** |
|:------------ |
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://gridap.github.io/Gridap.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gridap.github.io/Gridap.jl/dev) |
|**Build Status** |
| [![Build Status](https://travis-ci.com/gridap/Gridap.jl.svg?branch=master)](https://travis-ci.com/gridap/Gridap.jl) [![Codecov](https://codecov.io/gh/gridap/Gridap.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gridap/Gridap.jl) |
| **Community** |
| [![Join the chat at https://gitter.im/Gridap-jl/community](https://badges.gitter.im/Gridap-jl/community.svg)](https://gitter.im/Gridap-jl/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) |
| **Citation** |
| [![DOI](https://joss.theoj.org/papers/10.21105/joss.02520/status.svg)](https://doi.org/10.21105/joss.02520) |



Gridap provides a set of tools for the grid-based approximation of partial differential equations (PDEs) written in the
[Julia programming language](https://julialang.org/). The main motivation behind the development of this library is to provide an easy-to-use framework for the development of complex PDE solvers in a dynamically typed style without sacrificing the performance of statically typed languages. The library currently supports linear and nonlinear PDE systems for scalar and vector fields, single and multi-field problems, conforming and nonconforming finite element discretizations, on structured and unstructured meshes of simplices and hexahedra.

## Documentation

- [**STABLE**](https://gridap.github.io/Gridap.jl/stable) &mdash; **Documentation for the most recently tagged version of Gridap.jl.**
- [**DEVEL**](https://gridap.github.io/Gridap.jl/dev) &mdash; *Documentation for the in-development version of Gridap.*

## Tutorials

A hands-on user-guide to the library is available as a set of [tutorials](https://github.com/gridap/Tutorials). They are available as Jupyter notebooks and html pages.

## Installation

Gridap is a registered package in the official [Julia package registry](https://github.com/JuliaRegistries/General).  Thus, the installation of Gridap is straight forward using the [Julia's package manager](https://julialang.github.io/Pkg.jl/v1/). Open the Julia REPL, type `]` to enter package mode, and install as follows
```julia
pkg> add Gridap
```

## Plugins

- [GridapGmsh](https://github.com/gridap/GridapGmsh.jl) Generate a FE mesh with [GMSH](www.gmsh.info) and use it in Gridap.
- [GridapMakie](https://github.com/gridap/GridapMakie.jl) Makie plotting recipies for Gridap.
- [GridapPardiso](https://github.com/gridap/GridapPardiso.jl) Use the [Intel Pardiso MKL direct sparse solver](https://software.intel.com/en-us/mkl-developer-reference-fortran-intel-mkl-pardiso-parallel-direct-sparse-solver-interface) in Gridap.
- [GridapEmbedded](https://github.com/gridap/GridapEmbedded.jl) Embedded finite elements in Julia.
- [GridapODEs](https://github.com/gridap/GridapODEs.jl) Gridap support for time-dependent PDEs.
- [GridapDistributed](https://github.com/gridap/GridapDistributed.jl) Distributed-memory extension of Gridap.

## Examples

These are some popular PDEs solved with the Gridap library. Examples taken from the [Gridap Tutorials](https://github.com/gridap/Tutorials).

| ![](https://gridap.github.io/Tutorials/dev/assets/poisson/fig_uh.png)   |  ![](https://gridap.github.io/Tutorials/dev/assets/elasticity/disp_ux_40.png) | ![](https://gridap.github.io/Tutorials/dev/assets/hyperelasticity/neo_hook_3d.png)  | ![](https://gridap.github.io/Tutorials/dev/assets/p_laplacian/sol-plap.png)  |
|:-------------:|:-------------:|:-----:|:----:|
| [Poisson equation](https://gridap.github.io/Tutorials/dev/pages/t001_poisson/) |  [Linear elasticity](https://gridap.github.io/Tutorials/dev/pages/t003_elasticity/) |  [Hyper-elasticity](https://gridap.github.io/Tutorials/dev/pages/t005_hyperelasticity/)  | [p-Laplacian](https://gridap.github.io/Tutorials/dev/pages/t004_p_laplacian/)   |
| ![](https://gridap.github.io/Tutorials/dev/assets/dg_discretization/jump_u.png) | ![](https://gridap.github.io/Tutorials/dev/assets/darcy/darcy_results.png) |![](https://gridap.github.io/Tutorials/dev/assets/inc_navier_stokes/ins_solution.png) | ![](https://gridap.github.io/Tutorials/dev/assets/isotropic_damage/damage_end.png) |
| [Poisson eq. with DG](https://gridap.github.io/Tutorials/dev/pages/t006_dg_discretization/)  |  [Darcy eq. with RT](https://gridap.github.io/Tutorials/dev/pages/t007_darcy/)  |  [Incompressible Navier-Stokes](https://gridap.github.io/Tutorials/dev/pages/t008_inc_navier_stokes/)  | [Isotropic damage](https://gridap.github.io/Tutorials/dev/pages/t010_isotropic_damage/)  |

 ## Gridap community

Join to our [gitter](https://gitter.im/Gridap-jl/community) chat to ask questions and interact with the Gridap community.

## Contributing to Gridap

Gridap is a collaborative project open to contributions. If you want to contribute, please take into account:

  - Before opening a PR with a significant contribution, contact the project administrators, e.g., by writing a message in [our gitter chat](https://gitter.im/Gridap-jl/community) or by opening an issue describing what you are willing to implement. Wait for feed-back.
  - Carefully read and follow the instructions in the [CONTRIBUTING.md](https://github.com/gridap/Gridap.jl/blob/master/CONTRIBUTING.md) file.
  - Carefully read and follow the instructions in the [CODE_OF_CONDUCT.md](https://github.com/gridap/Gridap.jl/blob/master/CODE_OF_CONDUCT.md) file.
  - Open a PR with your contribution.

Want to help? We have a number of [issues waiting for help](https://github.com/gridap/Gridap.jl/labels/help%20wanted). You can start contributing to the Gridap project by solving some of those issues.


## How to cite Gridap

In order to give credit to the `Gridap` contributors, we simply ask you to cite the refence below in any publication in which you have made use of `Gridap` packages:

```
@article{Badia2020,
  doi = {10.21105/joss.02520},
  url = {https://doi.org/10.21105/joss.02520},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2520},
  author = {Santiago Badia and Francesc Verdugo},
  title = {Gridap: An extensible Finite Element toolbox in Julia},
  journal = {Journal of Open Source Software}
}
```

## Contact


Please, contact the project administrators, [Santiago Badia](mailto:santiago.badia@monash.edu) and [Francesc Verdugo](mailto:fverdugo@cimne.upc.edu), for further questions about licenses and terms of use.

