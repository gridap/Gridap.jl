
```@meta
CurrentModule = Gridap.Geometry
```
# Gridap.Geometry

```@docs
Geometry
```

## Triangulations

### Interface

```@docs
Triangulation
get_cell_coordinates(trian::Triangulation)
get_reffes(trian::Triangulation)
get_cell_type(trian::Triangulation)
test_triangulation
num_cells(trian::Triangulation)
num_cell_dims(::Triangulation{Dc,Dp}) where {Dc,Dp}
num_point_dims(::Triangulation{Dc,Dp}) where {Dc,Dp}
num_dims(g::Triangulation{Dc}) where Dc
is_affine(trian::Triangulation)
has_straight_faces(trian::Triangulation)
get_cell_reffes(trian::Triangulation)
get_cell_shapefuns(trian::Triangulation)
get_cell_map(trian::Triangulation)
```
## ConformingTriangulations

### Interface

```@docs
ConformingTriangulation
get_node_coordinates(trian::ConformingTriangulation)
get_cell_nodes(trian::ConformingTriangulation)
test_conforming_triangulation
num_nodes(trian::ConformingTriangulation)
```

### UnstructuredGrids

```@docs
UnstructuredGrid
UnstructuredGrid(
  node_coordinates::Vector{Point{Dp,Tp}},
  cell_nodes::Table{Ti},
  reffes::Vector{<:NodalReferenceFE{Dc}},
  cell_types::Vector) where {Dc,Dp,Tp,Ti}
UnstructuredGrid(trian::ConformingTriangulation)
UnstructuredGrid(reffe::NodalReferenceFE)
```

### CartesianGrids

```@docs
CartesianGrid
CartesianGrid(desc::CartesianDescriptor{D,T,F}) where {D,T,F}
CartesianGrid(domain,partition,map::Function)
get_cartesian_descriptor(a::CartesianGrid)
CartesianDescriptor
CartesianDescriptor(origin,sizes,partition,map::Function)
CartesianDescriptor(domain,partition,map::Function=identity)
```

## DiscreteModels

### Interface

```@docs
DiscreteModel
get_faces(g::DiscreteModel,dimfrom::Integer,dimto::Integer)
get_vertex_node(g::DiscreteModel)
get_node_face_owner(g::DiscreteModel)
get_face_nodes(g::DiscreteModel,d::Integer)
get_isboundary_face(g::DiscreteModel,d::Integer)
get_face_reffe_type(g::DiscreteModel,d::Integer)
get_face_polytope_type(g::DiscreteModel,d::Integer)
get_reffes(::Type{<:ReferenceFE{d}},g::DiscreteModel) where d
get_polytopes(::Type{<:Polytope{d}},g::DiscreteModel) where d
num_faces(g::DiscreteModel,d::Integer)
num_facets(g::DiscreteModel)
num_edges(g::DiscreteModel)
num_vertices(g::DiscreteModel)
get_vertex_coordinates(g::DiscreteModel)
get_dimranges(g::DiscreteModel)
get_offsets(g::DiscreteModel)
get_offset(g::DiscreteModel,d::Integer)
get_facedims(g::DiscreteModel)
get_cell_to_faces(g::DiscreteModel)
get_isboundary_face(g::DiscreteModel)
get_isboundary_node(g::DiscreteModel)
test_discrete_model
```


### UnstructuredDiscreteModel

!!! note
    Work-in-progress

