# Gridap.Adaptivity

```@meta
CurrentModule = Gridap.Adaptivity
```

The adaptivity module provides a framework to work with adapted (refined/coarsened/mixed) meshes.

It provides

- A generic interface to represent adapted meshes and a set of tools to work with Finite Element spaces defined on them. In particular, moving `CellFields` between parent and child meshes.
- Particular implementations for conformally refining/coarsening 2D/3D meshes using several well-known strategies. In particular, Red-Green refinement and longest-edge bisection.

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

The module provides a `refine` method for `UnstructuredDiscreteModel`. The method takes a string `refinement_method`
that determines the refinement strategy to be used. The following strategies are available:

- `"red_green"` :: Red-Green refinement, default.
- `"nvb"` :: Longest-edge bisection (only for meshes of TRIangles)
- `"barycentric"` :: Barycentric refinement (only for meshes of TRIangles)
- `"simplexify"` :: Simplexify refinement. Same resulting mesh as the `simplexify` method, but keeps track of the parent-child relationships.

Additionally, the method takes a kwarg `cells_to_refine` that determines which cells will be refined.
Possible input types are:

- `Nothing` :: All cells get refined.
- `AbstractArray{<:Bool}` of size `num_cells(model)` :: Only cells such that `cells_to_refine[iC] == true` get refined.
- `AbstractArray{<:Integer}` :: Cells for which `gid ∈ cells_to_refine` get refined

The algorithms try to respect the `cells_to_refine` input as much as possible, but some additional cells might get refined in order to guarantee that the mesh remains conforming.

```julia
  function refine(model::UnstructuredDiscreteModel;refinement_method="red_green",kwargs...)
    [...]
  end
```

## CartesianDiscreteModel refining

The module provides a `refine` method for `CartesianDiscreteModel`. The method takes a `Tuple` of size `Dc` (the dimension of the model cells) that will determine how many times cells will be refined in each direction. For example, for a 2D model, `refine(model,(2,3))` will refine each QUAD cell into a 2x3 grid of cells.

```julia
  function refine(model::CartesianDiscreteModel{Dc}, cell_partition::Tuple) where Dc
    [...]
  end
```

## Macro Finite-Elements

The module also provides support for macro finite-elements. From an abstract point of view, a macro finite-element is a finite-element defined on a refined polytope, where polynomial basis are defined on each of the subcells (creating a broken piece-wise polynomial space on the original polytope). From Gridap's point of view, a macro finite-element is a `ReferenceFE` defined on a `RefinementRule` from an array of `ReferenceFE`s defined on the subcells.

Although there are countless combinations, here are two possible applications:

- Linearized High-Order Lagrangian FESpaces: These are spaces which have the same DoFs as a high-order Lagrangian space, but where the basis functions are linear on each subcell.
- Barycentrically-refined elements for Stokes-like problems: These are spaces where the basis functions for velocity are defined on the barycentrically-refined mesh, whereas the basis functions for pressure are defined on the original cells. This allows for exact so-called Stokes sequences (see [here](https://arxiv.org/abs/2002.02051)).

The API is given by the following methods:

```@docs
  MacroReferenceFE
```

## Notes for users

Most of the tools provided by this module are showcased in the tests of the module itself, as well as the following tutorial (coming soon).

However, we want to stress a couple of key performance-critical points:

- The refining/coarsening routines are not optimized for performance. In particular, they are not parallelized.
  If you require an optimized/parallel implementation, please consider leveraging specialised meshing libraries. For instance, we provide an implementation of `refine/coarsen` using P4est in the [GridapP4est.jl](https://github.com/gridap/GridapP4est.jl) library.

- Although the toolbox allows you to evaluate `CellFields` defined on both fine/coarse meshes on their parent/children mesh, both directions of evaluation are not equivalent. As a user, you should always try to evaluate/integrate on the finest mesh for maximal performance. Evaluating a fine `CellField` on a coarse mesh relies on local tree searches, and is therefore a very expensive operation that should be avoided whenever possible.

## Notes for developers

### RefinementRule API

Given a `RefinementRule`, the library provides a set of methods to compute the mappings between parent (coarse) face ids and child (fine) face ids (and vice-versa).

The most basic information (that can directly be hardcoded in the RefinementRule for performance) are the mappings between parent face ids and child face ids. These are provided by:

```@docs
get_d_to_face_to_child_faces
get_d_to_face_to_parent_face
```

On top of these two basic mappings, a whole plethora of additional topological mappings can be computed.
These first set of routines extend the ReferenceFEs API to provide information on the face-to-node mappings and permutations:

```@docs
ReferenceFEs.get_face_vertices
ReferenceFEs.get_face_coordinates
ReferenceFEs.get_vertex_permutations
ReferenceFEs.get_face_vertex_permutations
```

We also provide face-to-face maps:

```@docs
get_cface_to_num_own_ffaces
get_cface_to_own_ffaces
get_cface_to_ffaces
get_cface_to_own_ffaces_to_lnodes
get_cface_to_ffaces_to_lnodes
get_cface_to_fface_permutations
aggregate_cface_to_own_fface_data
get_face_subface_ldof_to_cell_ldof
```

### AdaptivityGlue API

```@docs
get_n2o_reference_coordinate_map
get_old_cell_refinement_rules
get_new_cell_refinement_rules
get_d_to_fface_to_cface
n2o_reindex
o2n_reindex
```

### New-to-old field evaluations

When a cell is refined, we need to be able to evaluate the fields defined on the children cells on the parent cell. To do so, we bundle the fields defined on the children cells into a new type of `Field` called `FineToCoarseField`. When evaluated on a `Point`, a `FineToCoarseField` will select the child cell that contains the `Point` and evaluate the mapped point on the corresponding child field.

```@docs
FineToCoarseField
FineToCoarseDofBasis
FineToCoarseRefFE
```
