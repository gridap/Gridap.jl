

```@meta
CurrentModule = Gridap.Adaptivity
```

# Gridap.Adaptivity

The adaptivity module provides a framework to work with adapted (refined/coarsened/mixed) meshes.

It provides

  1) A generic interface to represent adapted meshes and a set of tools to work with Finite Element spaces defined on them. In particular, moving `CellFields` between different levels of the hierarchy.
  2) Particular implementations for conformally refining/coarsening 2D/3D meshes using several well-known strategies. In particular, Red-Green refinement and longest-edge bisection (TRI only).

## Interface

The following types are defined in the module:

```@docs
RefinementRule
AdaptivityGlue
AdaptedDiscreteModel
AdaptedTriangulation
```

The high-level interface is provided by the following methods:

```@docs
refine
coarsen
adapt
```

## Edge-Based refinement

Provides a `refine` method for `UnstructuredDiscreteModel`. The method takes a string `refinement_method`
that determines the refinement startegy to be used. The following strategies are available:

  - `"red_green"` :: Red-Green refinement, default.
  - `"nvb"` :: Longest-edge bisection (only for meshes of TRIangles)

Additionally, the method takes a kwarg `cells_to_refine` that determines which cells will be refined. 
Possible input types are:

  - `Nothing` :: All cells get refined.
  - `AbstractArray{<:Bool}` of size `num_cells(model)` :: Only cells such that `cells_to_refine[iC] == true` get refined.
  - `AbstractArray{<:Integer}` :: Cells for which `gid ∈ cells_to_refine` get refined

The algorithms try to respect the `cells_to_refine` input as much as possible, but some additional cells
might get refined in order to guarantee that the mesh remains conforming.

```julia
  function refine(model::UnstructuredDiscreteModel;refinement_method="red_green",kwargs...)
    [...]
  end
```

## CartesianDiscreteModel refining

Provides a `refine` method for `CartesianDiscreteModel`. The method takes a `Tuple` of size `Dc`
(the dimension of the model cells) that will determine how many times cells will be refined in
each direction. For example, for a 2D model, `refine(model,(2,3))` will refine each QUAD cell into
a 2x3 grid of cells.

```julia
  function refine(model::CartesianDiscreteModel{Dc}, cell_partition::Tuple) where Dc
    [...]
  end
```

## Notes for users

Most of the tools provided by this module are showcased in the tests of the module itself, as well as the following tutorial (comming soon).

However, we want to stress a couple of key performance-critical points:

- The refining/coarsening routines are not optimized for performance. In particular, they are not parallelized.
  If you require an optimized/parallel implementation, please consider leveraging spetialised meshing libraries. For instance, we provide an implementation of `refine/coarsen` using P4est in the [GridapP4est.jl](https://github.com/gridap/GridapP4est.jl) library.

- Although the toolbox allows you to evaluate `CellFields` defined on both fine/coarse meshes on their parent/children mesh, both directions of evaluation are not equivalent. As a user, you should always try to evaluate/integrate on the finest mesh for maximal performance. Evaluating a fine `CellField` on a coarse mesh relies on local tree searches, and is therefore a very expensive operation that should be avoided whenever possible.
