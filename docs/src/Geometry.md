
```@meta
CurrentModule = Gridap.Geometry
```

# Gridap.Geometry

```@contents
Pages = ["Geometry.md"]
Depth = 3
```

## DiscreteModels

In Gridap, a `DiscreteModel` is the main object representing a discretized domain, where the finite element problem is defined adn solved. It has three main components: 

- A `GridTopology`, which defines the topology of the mesh. That is the polytopes that make up the mesh as well as their connectivity. This object also provides the connectivity between all d-dimensional entities of the mesh, such as vertices, edges, faces, and cells.
- A `Grid`, which defines the geometric properties of the mesh. That is, the coordinates of the vertices and the geometric map between the reference space of each polytope and the physical space.
- A `FaceLabeling`, which classifies all the d-dimensional entities of the mesh into (non-overlapping) physical entities. These entities are then associated one or multiple tags which can be used to impose boundary conditions or define physical regions on the mesh.

To get a deeper understanding of these three components, we encourage the reader to check the low-level [tutorial on geometry](https://gridap.github.io/Tutorials/stable/). The API for each of these components is documented in the following sections.

### GridTopology API

```@autodocs
Modules = [Geometry,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/GridTopologies.jl","/GridTopologyMocks.jl"]
```

### Grid API

```@autodocs
Modules = [Geometry,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/Grids.jl","/GridMocks.jl"]
```

### FaceLabeling API

```@autodocs
Modules = [Geometry,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/FaceLabelings.jl"]
```

### DiscreteModel API

```@autodocs
Modules = [Geometry,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/DiscreteModels.jl","/DiscreteModelMocks.jl"]
```

## Types of models

### UnstructuredDiscreteModels

These are the main types of models used for unstructured meshes. They have their own type of topology and grid.

```@autodocs
Modules = [Geometry,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/UnstructuredDiscreteModels.jl","/UnstructuredGrids.jl","/UnstructuredGridTopologies.jl"]
```

### CartesianDiscreteModels

Cartesian models are specific structures to represent cartesian domains, possibly mapped by a prescribed function. They offer optimizations with respect to the unstructured models and are quite lighweight memory-wise. They implement their own type of grid and can be converted to unstructured models if needed.

```@autodocs
Modules = [Geometry,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/CartesianDiscreteModels.jl","/CartesianGrids.jl"]
```

### PolytopalDiscreteModels

Polytopal models are used to represent meshes made of arbitrarily shaped polytopes. They are directly build on the physical domain and therefore do require a geometric map. They can only be used with polytopal methods such as DG, HDG or HHO. More expensive than regular unstructured models, but also more flexible. They implement their own type of grid and topology.

```@autodocs
Modules = [Geometry,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/PolytopalDiscreteModels.jl"]
```

### Other models

```@autodocs
Modules = [Geometry,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/DiscreteModelPortions.jl","/MappedDiscreteModels.jl","/GridPortions.jl"]
```

## Triangulations

Given a model, which holds information for all dimensions and on the whole domain, we can take d-dimensional slices of it (or part of it) where we will define our finite element spaces and/or integrate our weakforms. These slices are called `Triangulations`.

### Body-Fitted triangulations

The most basic type of triangulation is the `BodyFittedTriangulation`. It represents a `d`-dimensional set of faces attached to `d`-dimensional faces of the model. In particular, it is the type used to represent bulk triangulations (which contain a subset of the cells of the model).

```@autodocs
Modules = [Geometry,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/Triangulations.jl"]
```

### Boundary and Skeleton triangulations

To perform integration of bulk variables on a set of faces (for instance the integration of Neumann terms), we somehow need to link these faces to one of their neighboring cells (otherwise the bulk variables are not uniquely valued on the faces).
This is done through the `BoundaryTriangulation` object. It is generally used for integrals on the physical boundary (where faces have a single neighboring cell), but offer quite a lot more flexibility. Conceptually, a `BoundaryTriangulation` is given by a `BodyFittedTriangulation` of faces and a `FaceToCellGlue` that links each face to the selected neighboring cell.
To perform integration on interior faces, which have two neighboring cells, we can use the `SkeletonTriangulation` object. Each `SkeletonTriangulation` is a wrapper around two `BoundaryTriangulations`, one for each side of the face. It provides straighforward ways to performs jumps and means over the interior faces.

```@autodocs
Modules = [Geometry,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/BoundaryTriangulations.jl","/SkeletonTriangulations.jl","/CompressedCellArrays.jl"]
```

### Patch triangulations

We provide an API to integrate and solve problems on patches of cells. This API revolves around two objects: The `PatchTopology` defines the topology of the patches, and holds information on all the d-dimensional entities belonging to each patch. One can then take d-dimensional slices of these patches using `PatchTriangulations`.

```@autodocs
Modules = [Geometry,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/PatchTriangulations.jl"]
```

### Other triangulations

```@autodocs
Modules = [Geometry,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/AppendedTriangulations.jl"]
```
