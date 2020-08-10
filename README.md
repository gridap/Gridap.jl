# Gridap

| **Documentation** |
|:------------ |
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://gridap.github.io/Gridap.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gridap.github.io/Gridap.jl/dev) [![Arxiv: 1910.01412](https://img.shields.io/badge/cs.MS-arXiv%3A1910.01412-B31B1B.svg)](https://arxiv.org/abs/1910.01412) |
|**Build Status** |
| [![Build Status](https://travis-ci.com/gridap/Gridap.jl.svg?branch=master)](https://travis-ci.com/gridap/Gridap.jl) [![Codecov](https://codecov.io/gh/gridap/Gridap.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gridap/Gridap.jl) |
| **Community** |
| [![Join the chat at https://gitter.im/Gridap-jl/community](https://badges.gitter.im/Gridap-jl/community.svg)](https://gitter.im/Gridap-jl/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) |
| **Citation** |
| [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3934468.svg)](https://doi.org/10.5281/zenodo.3934468) |



Gridap provides a set of tools for the grid-based approximation of partial differential equations (PDEs) written in the
[Julia programming language](https://julialang.org/). The main motivation behind the development of this library is to provide an easy-to-use framework for the development of complex PDE solvers in a dynamically typed style without sacrificing the performance of statically typed languages. The library currently supports linear and nonlinear PDE systems for scalar and vector fields, single and multi-field problems, conforming and nonconforming finite element discretizations, on structured and unstructured meshes of simplices and hexahedra.

## Documentation

- [**STABLE**](https://gridap.github.io/Gridap.jl/stable) &mdash; **Documentation for the most recently tagged version of Gridap.jl.**
- [**DEVEL**](https://gridap.github.io/Gridap.jl/dev) &mdash; *Documentation for the in-development version of Gridap.*
- [**ARTICLE**](https://arxiv.org/abs/1910.01412) &mdash; F. Verdugo, S. Badia. A user-guide to Gridap -- grid-based approximation of partial differential equations in Julia. *arXiv*. 2019. [arXiv:1910.01412](https://arxiv.org/abs/1910.01412)

## Installation

Gridap is a registered package in the official [Julia package registry](https://github.com/JuliaRegistries/General).  Thus, the installation of Gridap is straight forward using the [Julia's package manager](https://julialang.github.io/Pkg.jl/v1/). Open the Julia REPL, type `]` to enter package mode, and install as follows
```julia
pkg> add Gridap
```

## Tutorials

A hands-on user-guide to the library is available as a set of [tutorials](https://github.com/gridap/Tutorials). They are available as Jupyter notebooks and html pages.

## Plugins

- [GridapGmsh](https://github.com/gridap/GridapGmsh.jl) Generate a FE mesh with [GMSH](www.gmsh.info) and use it in Gridap.
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


## How to cite Gridap

If you have used the Gridap library in a scientific publication, please cite the project as follows:

```
@software{gridap_project,
  author       = {Francesc Verdugo and
                  Santiago Badia and
                  Víctor Sande and
                  Alberto F. Martin and
                  Oriol Colomés and
                  Jesús Bonilla},
  title        = {Gridap.jl},
  year         = 2020,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.3934468},
  url          = {https://doi.org/10.5281/zenodo.3934468}
}
```

## Contact


Please, contact the project administrators, [Santiago Badia](mailto:santiago.badia@monash.edu) and [Francesc Verdugo](mailto:fverdugo@cimne.upc.edu), for further questions about licenses and terms of use.

