# Gridap

[![Build Status](https://travis-ci.com/gridap/Gridap.jl.svg?branch=master)](https://travis-ci.com/gridap/Gridap.jl)
[![Codecov](https://codecov.io/gh/gridap/Gridap.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gridap/Gridap.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://gridap.github.io/Gridap.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gridap.github.io/Gridap.jl/dev) [![Join the chat at https://gitter.im/Gridap-jl/community](https://badges.gitter.im/Gridap-jl/community.svg)](https://gitter.im/Gridap-jl/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


## What

Gridap provides a set of tools for the grid-based approximation of partial differential equations (PDEs) written in the
[Julia programming language](https://julialang.org/). The main motivation behind the development of this library is to provide an easy-to-use framework for the development of complex PDE solvers in a dynamically typed style without sacrificing the performance of statically typed languages. The library currently supports linear and nonlinear PDE systems for scalar and vector fields, single and multi-field problems, conforming and nonconforming finite element discretizations, on structured and unstructured meshes of simplices and hexahedra.

## How

For further info, visit the project documentation:

- [**STABLE**](https://gridap.github.io/Gridap.jl/stable) &mdash; **Documentation for the most recently tagged version of Gridap.jl.**
- [**DEVEL**](https://gridap.github.io/Gridap.jl/dev) &mdash; *Documentation for the in-development version of Gridap.jl.*

## Examples

These are some popular PDE systems solved with the Gridap library. Examples taken from the [Gridap Tutorials](https://github.com/gridap/Tutorials).

| ![](https://gridap.github.io/Tutorials/dev/assets/t001_poisson/fig_uh.png)   |  ![](https://gridap.github.io/Tutorials/dev/assets/t003_elasticity/disp_ux_40.png) | ![](https://gridap.github.io/Tutorials/dev/assets/t004_hyperelasticity/neo_hook_3d.png)  | ![](https://gridap.github.io/Tutorials/dev/assets/t0041_p_laplacian/sol-plap.png)  |
|:-------------:|:-------------:|:-----:|:----:|
| [Poisson equation](https://gridap.github.io/Tutorials/dev/pages/t001_poisson/) |  [Linear elasticity](https://gridap.github.io/Tutorials/dev/pages/t003_elasticity/) |  [Hyper-elasticity](https://gridap.github.io/Tutorials/dev/pages/t004_hyperelasticity/)  | [p-Laplacian](https://gridap.github.io/Tutorials/dev/pages/t0041_p_laplacian/)   |
| ![](https://gridap.github.io/Tutorials/dev/assets/t006_poisson_dg/jump_u.png) | ![](https://gridap.github.io/Tutorials/dev/assets/t007_darcy/darcy_results.png) |![](https://gridap.github.io/Tutorials/dev/assets/t008_inc_navier_stokes/ins_solution.png) | |
| [Poisson eq. with DG](https://gridap.github.io/Tutorials/dev/pages/t005_dg_discretization/)  |  [Darcy eq. with RT](https://gridap.github.io/Tutorials/dev/pages/t007_darcy/)  |  [Incompressible Navier-Stokes](https://gridap.github.io/Tutorials/dev/pages/t008_inc_navier_stokes/)  |  |
