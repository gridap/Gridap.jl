# Gridap.jl

Documentation of the Gridap library.

!!! note

     These documentation pages are under construction.

## Introduction

Gridap provides a set of tools for the grid-based approximation of
partial differential equations (PDEs) written in the
[Julia programming language](https://julialang.org/).
The main motivation behind the development of this library is to provide an easy-to-use framework for the development of complex PDE solvers in a dynamically typed style without sacrificing the performance of statically typed languages.
The library currently supports linear and nonlinear PDE systems for scalar and vector fields, single and multi-field problems, conforming and nonconforming finite element discretizations, on structured and unstructured meshes of simplices and hexahedra.

## How to use this documentation

* The first step for new users is to visit the [Getting Started](@ref) page.

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
  "Gridap.md",
  "Helpers.md",
  "Io.md",
  "Algebra.md",
  "Arrays.md",
  "TensorValues.md",
  "Fields.md",
  "Polynomials.md",
  "Integration.md",
  "ReferenceFEs.md",
  "Geometry.md",
  "CellData.md",
  "Visualization.md",
  "FESpaces.md",
  "MultiField.md",
  "ODEs.md",
  "Adaptivity.md",
  ]
```
