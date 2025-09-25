# Gridap.jl

Documentation of the Gridap library.

## Introduction

Gridap provides a set of tools for the grid-based approximation of partial differential equations (PDEs) written in the [Julia programming language](https://julialang.org/).
The main motivation behind the development of this library is to provide an easy-to-use framework for the development of complex PDE solvers in a dynamically typed style without sacrificing the performance of statically typed languages.
The library currently supports linear and nonlinear PDE systems for scalar and vector fields, single and multi-field problems, conforming and nonconforming finite element discretizations, on structured and unstructured meshes of simplices and hexahedra.

## How to use this documentation

* The first step for new users is to visit the [Getting Started](@ref) page.

* For a high-level overview of Gridap's capabilities, see [Gridap at a glance](@ref).

* To explore the extended capabilities available through companion packages, check out the [Gridap ecosystem](@ref).

* A set of tutorials written as Jupyter notebooks and html pages are available [here](https://github.com/gridap/Tutorials).

* The detailed documentation is in the [Manual](@ref) section.

* Guidelines for developers of the Gridap project is found in the [Gridap wiki](https://github.com/gridap/Gridap.jl/wiki) page.

## Julia educational resources

A basic knowledge of the Julia programming language is needed to use the Gridap package.
Here, one can find a list of resources to get started with this programming language.

* First steps to learn Julia form the [Gridap wiki](https://github.com/gridap/Gridap.jl/wiki/Start-learning-Julia) page.
* Official webpage [docs.julialang.org](https://docs.julialang.org/)
* Official list of learning resources [julialang.org/learning](https://julialang.org/learning/)

## Manual

```@contents
Pages = [
  "overview.md",
  "ecosystem.md",
  "modules/Helpers.md",
  "modules/Io.md",
  "modules/Algebra.md",
  "modules/Arrays.md",
  "modules/TensorValues.md",
  "modules/Fields.md",
  "modules/Polynomials.md",
  "Integration.md",
  "modules/ReferenceFEs.md",
  "modules/Geometry.md",
  "modules/CellData.md",
  "modules/Visualization.md",
  "modules/FESpaces.md",
  "modules/MultiField.md",
  "modules/ODEs.md",
  "modules/Adaptivity.md",
  ]
```
