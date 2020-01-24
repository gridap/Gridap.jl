
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
get_normal_vector(trian::Triangulation)
restrict(f::AbstractArray, trian::Triangulation)
reindex(f::AbstractArray,trian::Triangulation)
test_triangulation
num_cells(trian::Triangulation)
num_cell_dims(::Triangulation{Dc,Dp}) where {Dc,Dp}
num_point_dims(::Triangulation{Dc,Dp}) where {Dc,Dp}
num_dims(g::Triangulation{Dc}) where Dc
is_affine(trian::Triangulation)
is_first_order(trian::Triangulation)
get_cell_reffes(trian::Triangulation)
get_cell_shapefuns(trian::Triangulation)
get_cell_map(trian::Triangulation)
get_physical_coordinate(trian::Triangulation)
```

### TriangulationPortion

```@docs
TriangulationPortion
TriangulationPortion(oldtrian::Triangulation{Dc,Dp},cell_to_oldcell::Vector{Int}) where {Dc,Dp}
```

## BoundaryTriangulations

### Interface

```@docs
BoundaryTriangulation
BoundaryTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
BoundaryTriangulation(model::DiscreteModel,tags::Vector{Int})
get_volume_triangulation(trian::BoundaryTriangulation)
get_face_to_cell(trian::BoundaryTriangulation)
get_face_to_cell_map(trian::BoundaryTriangulation)
get_normal_vector(trian::BoundaryTriangulation)
test_boundary_triangulation
```

### GenericBoundaryTriangulations

```@docs
GenericBoundaryTriangulation
GenericBoundaryTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
```

## SkeletonTriangulations

```@docs
SkeletonTriangulation
SkeletonTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
get_volume_triangulation(trian::SkeletonTriangulation)
get_normal_vector(trian::SkeletonTriangulation)
```

### SkeletonPairs

```@docs
SkeletonPair
```

## Grids

### Interface

```@docs
Grid
OrientationStyle(a::Grid)
RegularityStyle(a::Grid)
get_node_coordinates(trian::Grid)
get_cell_nodes(trian::Grid)
test_grid
num_nodes(trian::Grid)
is_oriented(a::Grid)
is_regular(a::Grid)
Grid(reffe::LagrangianRefFE)
compute_linear_grid(reffe::LagrangianRefFE)
compute_reference_grid(reffe::LagrangianRefFE, nelems::Integer)
Grid(::Type{ReferenceFE{d}},p::Polytope) where d
GridTopology(grid::Grid)
```

### UnstructuredGrids

```@docs
UnstructuredGrid
UnstructuredGrid(
  node_coordinates::Vector{Point{Dp,Tp}},
  cell_nodes::Table{Ti},
  reffes::Vector{<:NodalReferenceFE{Dc}},
  cell_types::Vector,
  ::Val{B}=Val{false}()) where {Dc,Dp,Tp,Ti,B}
UnstructuredGrid(trian::Grid)
UnstructuredGrid(reffe::LagrangianRefFE)
UnstructuredGrid(::Type{ReferenceFE{d}},p::Polytope) where d
UnstructuredGrid(x::AbstractArray{<:Point})
UnstructuredGridTopology(grid::UnstructuredGrid)
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

### GridPortion

```@docs
GridPortion
GridPortion(oldgrid::Grid{Dc,Dp},cell_to_oldcell::Vector{Int}) where {Dc,Dp}
```

## FaceLabeling

```@docs
FaceLabeling
FaceLabeling(d_to_num_dfaces::Vector{Int})
num_dims(lab::FaceLabeling)
num_cell_dims(lab::FaceLabeling)
num_tags(lab::FaceLabeling)
num_entities(lab::FaceLabeling)
num_faces(lab::FaceLabeling,d::Integer)
num_faces(lab::FaceLabeling)
num_vertices(lab::FaceLabeling)
num_edges(lab::FaceLabeling)
num_facets(lab::FaceLabeling)
num_cells(lab::FaceLabeling)
get_face_entity(lab::FaceLabeling,d::Integer)
get_face_entity(lab::FaceLabeling)
get_tag_entities(lab::FaceLabeling,tag::Integer)
get_tag_entities(lab::FaceLabeling)
get_tag_name(lab::FaceLabeling,tag::Integer)
get_tag_name(lab::FaceLabeling)
get_tag_from_name(lab::FaceLabeling,name::String)
get_tag_from_name(lab::FaceLabeling)
get_tags_from_names(lab::FaceLabeling,names::Vector{String})
get_face_mask(labeling::FaceLabeling,tags::Vector{Int},d::Integer)
add_tag!(lab::FaceLabeling,name::String,entities::Vector{<:Integer})
add_tag_from_tags!(lab::FaceLabeling, name::String, tags::Vector{Int})
get_face_tag(labeling::FaceLabeling,tags,d::Integer)
get_face_tag_index(labeling::FaceLabeling,tags,d::Integer)
```

## GridTopology

### Interface

```@docs
GridTopology
OrientationStyle(a::GridTopology)
RegularityStyle(a::GridTopology)
get_faces(g::GridTopology,dimfrom::Integer,dimto::Integer)
get_polytopes(g::GridTopology)
get_cell_type(g::GridTopology)
get_vertex_coordinates(g::GridTopology)
test_grid_topology(top::GridTopology{Dc,Dp}) where {Dc,Dp}
num_cell_dims(::GridTopology{Dc,Dp}) where {Dc,Dp}
num_point_dims(::GridTopology{Dc,Dp}) where {Dc,Dp}
num_dims(g::GridTopology{Dc}) where Dc
num_faces(g::GridTopology,d::Integer)
num_cells(g::GridTopology)
num_facets(g::GridTopology)
num_edges(g::GridTopology)
num_vertices(g::GridTopology)
get_dimranges(g::GridTopology)
get_dimrange(g::GridTopology,d::Integer)
get_offsets(g::GridTopology)
get_offset(g::GridTopology,d::Integer)
get_facedims(g::GridTopology)
get_cell_faces(g::GridTopology)
compute_cell_faces(g::GridTopology)
get_face_vertices(g::GridTopology,d::Integer)
get_face_vertices(g::GridTopology)
compute_face_vertices(g::GridTopology)
get_cell_vertices(g::GridTopology)
is_simplex(p::GridTopology)
is_n_cube(p::GridTopology)
is_oriented(a::GridTopology)
is_regular(a::GridTopology)
get_reffaces(::Type{Polytope{d}}, g::GridTopology) where d
get_face_type(g::GridTopology,d::Integer)
compute_reffaces(::Type{Polytope{d}}, g::GridTopology) where d
get_reffaces(topo::GridTopology)
get_face_type(topo::GridTopology)
get_reffaces_offsets(topo::GridTopology)
compute_reffaces(g::GridTopology)
get_isboundary_face(g::GridTopology)
get_isboundary_face(g::GridTopology,d::Integer)
compute_isboundary_face(g::GridTopology)
compute_isboundary_face(g::GridTopology,d::Integer)
get_cell_permutations(top::GridTopology)
get_cell_permutations(top::GridTopology,d::Integer)
compute_cell_permutations(top::GridTopology)
compute_cell_permutations(top::GridTopology,d::Integer)
```
### UnstructuredGridTopology

```@docs
UnstructuredGridTopology
UnstructuredGridTopology(
  vertex_coordinates::Vector{<:Point},
  cell_vertices::Table,
  cell_type::Vector{<:Integer},
  polytopes::Vector{<:Polytope},
  orientation::Val{O}=Val{false}()) where O
UnstructuredGridTopology(
  vertex_coordinates::Vector{<:Point},
  d_to_dface_vertices::Vector{<:Table},
  cell_type::Vector{<:Integer},
  polytopes::Vector{<:Polytope},
  orientation::Val{O}=Val{false}()) where O
```

## DiscreteModels

### Interface

```@docs
DiscreteModel
get_grid(model::DiscreteModel)
get_grid_topology(model::DiscreteModel)
get_face_labeling(g::DiscreteModel)
test_discrete_model(model::DiscreteModel{Dc,Dp}) where {Dc,Dp}
num_dims(model::DiscreteModel)
num_cell_dims(model::DiscreteModel)
num_point_dims(model::DiscreteModel)
num_faces(g::DiscreteModel,d::Integer)
num_cells(g::DiscreteModel)
num_facets(g::DiscreteModel)
num_edges(g::DiscreteModel)
num_vertices(g::DiscreteModel)
num_nodes(g::DiscreteModel)
get_face_nodes(g::DiscreteModel,d::Integer)
get_face_nodes(g::DiscreteModel)
compute_face_nodes(model::DiscreteModel,d::Integer)
compute_face_nodes(model::DiscreteModel)
get_face_own_nodes(g::DiscreteModel,d::Integer)
get_face_own_nodes(g::DiscreteModel)
compute_face_own_nodes(model::DiscreteModel,d::Integer)
compute_face_own_nodes(model::DiscreteModel)
get_vertex_node(g::DiscreteModel)
compute_vertex_node(g::DiscreteModel)
get_node_face_owner(g::DiscreteModel)
compute_node_face_owner(g::DiscreteModel)
get_reffaces(::Type{ReferenceFE{d}},model::DiscreteModel) where d
get_face_type(g::DiscreteModel,d::Integer)
compute_reffaces(::Type{ReferenceFE{d}}, g::DiscreteModel) where d
get_reffaces(model::DiscreteModel)
get_face_type(model::DiscreteModel)
get_reffaces_offsets(model::DiscreteModel)
compute_reffaces(g::DiscreteModel)
Grid(::Type{ReferenceFE{d}},model::DiscreteModel) where d
Triangulation(::Type{ReferenceFE{d}},model::DiscreteModel) where d
get_triangulation(model::DiscreteModel)
get_polytopes(model::DiscreteModel)
```

### UnstructuredDiscreteModels

```@docs
UnstructuredDiscreteModel
UnstructuredDiscreteModel(trian::Grid)
```
### CartesianDiscreteModels

```@docs
CartesianDiscreteModel
CartesianDiscreteModel(desc::CartesianDescriptor{D,T,F}) where {D,T,F}
CartesianDiscreteModel(args...)
get_cartesian_descriptor(a::CartesianDiscreteModel)
```

## CellFields

### CellFieldLike interface
```@docs
CellFieldLike
get_array(cf::CellFieldLike)
get_cell_map(cf::CellFieldLike)
similar_object(cf::CellFieldLike,array::AbstractArray)
similar_object(cf1::CellFieldLike,cf2::CellFieldLike,array::AbstractArray)
gradient(cf::CellFieldLike)
grad2curl(cf::CellFieldLike)
test_cell_field_like
evaluate(cf::CellFieldLike,x)
length(cf::CellFieldLike)
```
### CellField interface

```@docs
CellField
test_cell_field
convert_to_cell_field(object::CellField,cell_map)
restrict(cf::CellField,trian::Triangulation)
```
### Concrete implementations

```@docs
GenericCellField
SkeletonCellField
get_cell_map(a::SkeletonCellField)
jump(sf::SkeletonCellField)
mean(sf::SkeletonCellField)
```
