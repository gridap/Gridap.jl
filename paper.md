---
title: 'Gridap: An extensible Finite Element toolbox in Julia'
tags:
  - julia
  - pdes
  - partial differential equations
  - finite elements
authors:
  - name: Santiago Badia
    orcid: 0000-0003-2391-4086
    affiliation: "1"
  - name: Francesc Verdugo
    orcid: 0000-0003-3667-443X
    affiliation: 2
affiliations:
 - name: School of Mathematics, Monash University, Clayton, Victoria, 3800, Australia.
   index: 1
 - name: Centre Internacional de Mètodes Numèrics en Enginyeria, Esteve Terrades 5, E-08860 Castelldefels, Spain.
   index: 2
date: 6 July 2020
bibliography: paper.bib
---

# Summary

Gridap is a new Finite Element (FE) framework, exclusively written in the Julia programming language, for the numerical simulation of a wide range of mathematical models governed by partial differential equations (PDEs). The library provides a feature-rich set of discretization techniques, including continuous and discontinuous FE methods with Lagrangian, Raviart-Thomas, or Nédélec interpolations, and supports a wide range of problem types including linear, nonlinear, single-field, and multi-field PDEs (see [@badia_fempar:_2017, Section 3] for a detailed presentation of the mathematical abstractions behind the implementation of these FE methods). Gridap is designed to help application experts to easily simulate real-world problems, to help researchers improve productivity when developing new FE-related techniques, and also for its usage in numerical PDE courses.

The main motivation behind Gridap is to find an improved balance between computational performance, user-experience, and work-flow productivity when working with FE libraries. Previous FE frameworks, e.g., FEniCS [@Alnaes2015] or Deal.II [@Bangerth2007] usually provides a high-level user front-end to facilitate the use of the library and a computational back-end to achieve performance. The user front-end is usually programmable in an interpreted language like Python, whereas the computational back-end is usually coded in a compiled language like C/C++ or Fortran. Users can benefit from the high-level front-end (i.e., for rapid prototyping) and simultaneously enjoy the performance of the compiled back-end. This approach reaches a compromise between performance and productivity when the back-end provides all the functionality required by the user. However, it does not satisfactorily address the needs of researchers on numerical methods willing to extend the library with new techniques or features. These extensions usually need to be done at the level of the computational back-end for performance reasons. Thus, the researcher is forced to develop a new code in a compiled language like C/C++ instead of benefiting from the productivity of scripting languages like Python, incurring serious productivity losses. In order to overcome this limitation, Gridap is fully implemented in the Julia programming language [@Bezanson2017]. Julia combines the performance of compiled languages with the productivity of interpreted ones by using type inference and just-in-time compilation to generate fast code. As a result, there is no need to use two different languages to write low-level performance code and high-level user interfaces. In addition, writing a FE library in Julia also allows one to leverage the feature-rich ecosystem of Julia libraries and exploit its excellent package manager. It permits a seamless coupling of Gridap with application-specific libraries, like optimization [@Dunning2017], an approximation of ordinary differential equations [@Rackauckas2017], or data science [@Innes2018], which can certainly boost the capabilities of a FE solver.

Another major feature of Gridap is that it is not a simple Julia translation of a standard object-oriented FE code. There are other FE libraries written in Julia that have been inspired by standard FE frameworks, see, e.g., JuAFEM [@Carlsson], whose interface resembles Deal.II.  In contrast,  Gridap adopts a novel software design that allows one to manipulate different types of data associated with the cells of the computational mesh in a convenient way. For instance, one can build an object representing the elemental matrices for all cells in the mesh using high-level API calls, without explicitly writing any for-loop. These objects representing data for all cells of the mesh are usually _lazy_, meaning that the underlying data is never stored for all cells in the mesh simultaneously. Instead, the value for a specific cell is computed on-the-fly when needed, which certainly reduces memory requirements.  This software design allows the library developers to hide assembly loops and other core computations from the user-code, leading to a very compact, user-friendly, syntax, while providing a high degree of flexibility for users to define their own weak forms.  A Poisson or Stokes problem can be solved with Gridap in 10-20 lines of code, as this example for the Poisson equation shows:

```julia
using Gridap
# Manufactured solutions
u(x) = x[1]^2 + x[2]
f(x) = -Δ(u)(x); g(x) = u(x)
# FE mesh (aka discrete model)
pmin = Point(0,0,0); pmax = Point(1,1,1)
cells=(8,8,8); order = 1
model = CartesianDiscreteModel(pmin, pmax, cells)
# FE Spaces
V0 = TestFESpace(model=model, reffe=:Lagrangian,
  valuetype=Float64, order=order,
  conformity=:H1, dirichlet_tags="boundary")
Ug = TrialFESpace(V0, g)
# Weak form
a(u,v) = ∇(u)⋅∇(v); l(v) = v*f
trian_Ω = Triangulation(model)
quad_Ω = CellQuadrature(trian_Ω, 2*order)
t_Ω = AffineFETerm(a,l,trian_Ω,quad_Ω)
# FE Problem and solution
op = AffineFEOperator(Ug,V0,t_Ω)
uh = solve(op)
# Output for visualization
writevtk(trian_Ω,"results",
  cellfields=["uh"=>uh,"grad_uh"=>∇(uh)])
```

Other FE packages like FEniCS also achieve such compact user interfaces, but in contrast to Gridap, they are based on a sophisticated compiler of variational forms [@Kirby2006], which generates, compiles and links a specialized C++ back-end for the problem at hand. One of the limitations of this approach is that the form compiler is a rigid system that is not designed to be extended by average users.

Gridap is an open-source project hosted at Github and distributed under an MIT license. The source code for Gridap has been archived to Zenodo with the linked DOI [10.5281/zenodo.3934468](https://doi.org/10.5281/zenodo.3934468).

# References

